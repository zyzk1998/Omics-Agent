"""
组学下游算子共享：可执行文件嗅探、参考序列解析、subprocess 日志裁剪与缺依赖时的诚实报错。

宿主具备 CLI 与输入资源时执行真实 subprocess；**禁止用仿真生物学指标静默顶替**。
重型全长管线仍推荐 Worker/TaaS。
"""
from __future__ import annotations

import glob
import logging
import os
import shutil
import subprocess
import tempfile
import traceback
from typing import Any, Callable, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# 与 Runner 微缩参考兜底一致（仅用于诊断展示，避免循环 import）
_MINIMAL_HG38_FALLBACK_PATH = "/tmp/omics_test_data/ref/chr21.fa"
_CONTAINER_DEFAULT_OMICS_REF = "/app/references"


def exe_available(name: str) -> bool:
    return resolve_cli_exe(name) is not None


def _is_executable_file(path: str) -> bool:
    try:
        return os.path.isfile(path) and os.access(path, os.X_OK)
    except OSError:
        return False


def _glob_conda_cli_candidates(base_exe: str) -> List[str]:
    """跨 Conda 安装布局的模糊 bin 扫描（不依赖进程 PATH）。"""
    patterns = (
        "/home/*/*conda*/envs/*/bin",
        "/opt/*conda*/envs/*/bin",
        "/root/*conda*/envs/*/bin",
    )
    found: List[str] = []
    for base_dir in patterns:
        for match in glob.glob(os.path.join(base_dir, base_exe)):
            if _is_executable_file(match):
                found.append(match)
    return sorted(set(found))


def _explicit_conda_cli_candidates(base_exe: str) -> List[str]:
    """常见容器 / 宿主机固定前缀（绝对路径穿透，绕过受限 PATH）。"""
    home = os.path.expanduser("~")
    env_names = ("omics-real", "base")
    roots: List[str] = []
    for env in env_names:
        roots.extend(
            [
                f"/opt/omics-conda/envs/{env}/bin",
                f"/opt/omics-conda/bin",
                f"/opt/conda/envs/{env}/bin",
                f"/opt/conda/bin",
                os.path.join(home, "miniconda3", "envs", env, "bin"),
                os.path.join(home, "miniforge3", "envs", env, "bin"),
                os.path.join(home, "anaconda3", "envs", env, "bin"),
                f"/root/miniconda3/envs/{env}/bin",
                f"/root/miniforge3/envs/{env}/bin",
            ]
        )
    roots.extend(
        [
            "/usr/local/bin",
            "/usr/bin",
            "/bin",
            "/opt/bin",
            "/opt/homebrew/bin",
            os.path.join(home, "miniconda3", "bin"),
            os.path.join(home, "miniforge3", "bin"),
            os.path.join(home, "anaconda3", "bin"),
            os.path.join(home, "micromamba", "bin"),
        ]
    )
    return [os.path.join(r, base_exe) for r in roots]


def resolve_cli_exe(exe_name: str) -> Optional[str]:
    """
    解析生信 CLI：优先 PATH（shutil.which），其次固定 Conda/系统前缀与 glob 深度扫描。
    命中后返回**绝对路径**，供 subprocess 直接调用，绕过 API 进程未继承 shell PATH 的限制。
    """
    name = (exe_name or "").strip()
    if not name:
        return None
    base_exe = os.path.basename(name)

    def _first_executable(candidates: List[str]) -> Optional[str]:
        for cand in candidates:
            if _is_executable_file(cand):
                return os.path.abspath(cand)
        return None

    w = shutil.which(base_exe)
    if w and _is_executable_file(w):
        return os.path.abspath(w)

    hit = _first_executable(_explicit_conda_cli_candidates(base_exe))
    if hit:
        return hit

    glob_hits = _glob_conda_cli_candidates(base_exe)
    if glob_hits:
        return glob_hits[0]

    search_roots: List[str] = []
    for env_root_name in ("miniconda3", "miniforge3", "anaconda3", "mambaforge", "micromamba"):
        envs = os.path.expanduser(f"~/{env_root_name}/envs")
        if not os.path.isdir(envs):
            continue
        try:
            sub = sorted(os.listdir(envs))
        except OSError:
            continue
        for env in sub[:48]:
            search_roots.append(os.path.join(envs, env, "bin"))

    return _first_executable([os.path.join(root, base_exe) for root in search_roots])


def clip_omics_log(text: str, max_len: int = 16000) -> str:
    t = (text or "").strip()
    if len(t) <= max_len:
        return t
    return t[: max_len - 48] + "\n\n… [log truncated] …\n"


def omics_frontend_error_from_exception(ctx: str, exc: BaseException) -> Dict[str, Any]:
    """不可预期异常：原样 traceback 交给前端 Markdown。"""
    from .omics_mock_ui import IMG_DNA_HELIX, attach_visual_contract, simple_rows_table

    tb = traceback.format_exc()
    md = (
        f"### 组学工具执行异常（{ctx}）\n\n"
        f"- **类型**: `{type(exc).__name__}`\n"
        f"- **消息**: `{exc}`\n\n"
        "<details><summary>完整 traceback（服务端捕获）</summary>\n\n```\n"
        f"{clip_omics_log(tb, 12000)}"
        "\n```\n</details>\n"
    )
    out: Dict[str, Any] = {
        "status": "error",
        "message": f"{ctx}: {exc}",
        "output_path": "",
        "file_path": "",
    }
    return attach_visual_contract(
        out,
        markdown=md,
        image_urls=[IMG_DNA_HELIX],
        table_data=simple_rows_table(
            ("field", "value"),
            [
                {"field": "error_type", "value": type(exc).__name__},
                {"field": "context", "value": ctx[:240]},
            ],
        ),
    )


def omics_subprocess_failed(
    ctx: str,
    cmd: List[str],
    proc: subprocess.CompletedProcess,
    *,
    tool_id: str = "",
) -> Dict[str, Any]:
    """CLI 非零退出：stdout/stderr 写入 Markdown。"""
    from .omics_mock_ui import IMG_DNA_HELIX, attach_visual_contract, simple_rows_table

    stdout = getattr(proc, "stdout", None)
    stderr = getattr(proc, "stderr", None)
    sout = (
        stdout.decode("utf-8", "replace")
        if isinstance(stdout, bytes)
        else (stdout or "")
    )
    serr = (
        stderr.decode("utf-8", "replace")
        if isinstance(stderr, bytes)
        else (stderr or "")
    )
    from .omics_genomics_report_ui import assemble_genomics_step_markdown

    cmd_s = " ".join(str(c) for c in cmd[:48])
    if len(cmd) > 48:
        cmd_s += " …"
    fake_cp = subprocess.CompletedProcess(cmd, proc.returncode, sout, serr)
    md = assemble_genomics_step_markdown(
        title=f"CLI 非零退出（{ctx}）",
        metrics=[
            ("returncode", str(proc.returncode)),
            ("stderr 长度", f"{len(serr)} 字符"),
            ("stdout 长度", f"{len(sout)} 字符"),
        ],
        body_md=f"- **命令**: `{cmd_s}`",
        cmd=cmd,
        cp=fake_cp,
    )
    out: Dict[str, Any] = {
        "status": "error",
        "message": f"{ctx} subprocess exit {proc.returncode}",
        "output_path": "",
        "file_path": "",
    }
    if tool_id:
        out["tool_id"] = tool_id
    return attach_visual_contract(
        out,
        markdown=md,
        image_urls=[IMG_DNA_HELIX],
        table_data=simple_rows_table(
            ("stream", "length_chars"),
            [
                {"stream": "stderr", "length_chars": str(len(serr))},
                {"stream": "stdout", "length_chars": str(len(sout))},
            ],
        ),
        tool_id=tool_id or None,
    )


def omics_host_prerequisite_blocked(
    *,
    context: str,
    tool_id: str,
    title: str,
    checks: List[Tuple[str, str]],
    extra_md: str = "",
) -> Dict[str, Any]:
    """缺依赖 / 缺输入：输出诊断表 + 完整 PATH / 参考路径，不伪造实验指标。"""
    from .omics_mock_ui import IMG_DNA_HELIX, attach_visual_contract, simple_rows_table

    rows_md = ""
    if checks:
        rows_md = "\n| 检查项 | 结果 |\n|--------|------|\n"
        for k, v in checks:
            rows_md += f"| {k} | {v} |\n"

    missing_tools = [
        k
        for k, v in checks
        if k and ("missing" in (v or "").lower() or (v or "").strip().lower() in ("blocked", "missing"))
    ]

    md = (
        f"### {title}\n\n"
        "**本步未输出实验测定结果。** 以下为宿主环境诊断（非生物学仿真）。\n\n"
        "请按仓库根目录 `scripts/install_omics_real_env.sh` 安装 CLI、参考序列与索引；"
        "并通过环境变量 `GIBH_REF_HG38`、`GIBH_BOWTIE2_HG38_INDEX` 等注入路径。\n"
        f"{rows_md}\n"
    )
    if extra_md.strip():
        md += "\n#### 补充说明\n\n" + extra_md.strip() + "\n"
    md += (
        _deep_env_dependency_markdown(
            title_line=f"工具 **`{tool_id}`**（上下文 `{context}`）前置检查未通过。",
            missing_tools=missing_tools,
        )
    )
    md += f"\n- **上下文**: `{context}`\n"

    out: Dict[str, Any] = {
        "status": "error",
        "message": title,
        "output_path": "",
        "file_path": "",
        "tool_id": tool_id,
        "fatal_env_prereq": True,
        "can_skip": False,
        "error_category": "host_env_prereq",
    }
    tbl_rows = [{"check": a, "status": b} for a, b in checks]
    if not tbl_rows:
        tbl_rows = [{"check": "diagnosis", "status": "blocked"}]
    return attach_visual_contract(
        out,
        markdown=md,
        image_urls=[IMG_DNA_HELIX],
        table_data=simple_rows_table(("check", "status"), tbl_rows),
        tool_id=tool_id,
    )


def run_omics_without_synthetic_fallback(
    real_try: Callable[[], Optional[Dict[str, Any]]],
    blocked_fn: Callable[[], Dict[str, Any]],
    ctx: str = "",
) -> Dict[str, Any]:
    """
    real_try 返回 dict 则视为终态（含 status:error）。
    None 表示调用 blocked_fn（诚实的环境说明）。
    """
    try:
        got = real_try()
        if got is not None:
            return got
    except Exception as exc:  # noqa: BLE001
        logger.exception("%s: real_try raised", ctx)
        return omics_frontend_error_from_exception(ctx or "omics", exc)
    return blocked_fn()


def write_temp_mock_artifact(prefix: str, message: str) -> str:
    fd, path = tempfile.mkstemp(prefix=f"{prefix}_", suffix=".mock.txt")
    os.close(fd)
    try:
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(message + "\n")
    except OSError as exc:
        logger.warning("mock artifact write failed: %s", exc)
    return path


def degraded_banner_md(missing: str) -> str:
    """遗留：仅用于仍调用旧文案的路径；新逻辑请用 omics_host_prerequisite_blocked。"""
    return (
        "> **运行环境提示**：未检测到 "
        + missing
        + "。\n\n"
    )


def get_omics_ref_dir() -> str:
    """统一参考库根目录：OMICS_REF_DIR > 容器默认 /app/references > 仓库 data/references。"""
    for cand in (
        os.environ.get("OMICS_REF_DIR", "").strip(),
        _CONTAINER_DEFAULT_OMICS_REF,
        os.path.join(os.getcwd(), "data", "references"),
    ):
        if cand and os.path.isdir(cand):
            return os.path.abspath(cand)
    return os.path.abspath(
        os.environ.get("OMICS_REF_DIR", "").strip() or _CONTAINER_DEFAULT_OMICS_REF
    )


def omics_integrated_infra_ready() -> bool:
    """data/references/manifest.json 存在表示基建脚本已跑通（参考 + 索引）。"""
    manifest = os.path.join(get_omics_ref_dir(), "manifest.json")
    return os.path.isfile(manifest)


def bwa_index_ready(fasta_path: str) -> bool:
    fa = (fasta_path or "").strip()
    return bool(fa and os.path.isfile(fa) and os.path.isfile(f"{fa}.bwt"))


def resolve_reference_fasta(reference_id: str) -> Optional[str]:
    """通过 reference_id、GIBH_REF_* 与 OMICS_REF_DIR 解析参考基因组 FASTA。"""
    rid = (reference_id or "").strip().lower().replace(" ", "")
    env_for_id = {
        "hg38": "GIBH_REF_HG38",
        "grch38": "GIBH_REF_HG38",
        "hg19": "GIBH_REF_HG19",
        "grch37": "GIBH_REF_HG19",
        "mm10": "GIBH_REF_MM10",
        "mm39": "GIBH_REF_MM39",
    }
    primary = env_for_id.get(rid)
    candidates: List[str] = []
    if primary:
        candidates.append(os.getenv(primary) or "")
    candidates.append(os.getenv("GIBH_REFERENCE_FASTA") or "")
    ref_root = get_omics_ref_dir()
    if rid in ("hg38", "grch38"):
        candidates.extend(
            [
                os.path.join(ref_root, "genomics", "hg38.fa"),
                os.path.join(ref_root, "genomics", "GRCh38.fa"),
            ]
        )
    elif rid in ("hg19", "grch37"):
        candidates.append(os.path.join(ref_root, "genomics", "hg19.fa"))
    for p in candidates:
        if p and os.path.isfile(p):
            return os.path.abspath(p)
    return None


def resolve_proteomics_search_fasta() -> Optional[str]:
    p = (os.getenv("GIBH_PROTEOMICS_FASTA") or "").strip()
    if p and os.path.isfile(p):
        return os.path.abspath(p)
    cand = os.path.join(get_omics_ref_dir(), "proteomics", "uniprot_mini.fasta")
    if os.path.isfile(cand):
        return os.path.abspath(cand)
    return None


def omics_infra_stage_success(
    *,
    tool_id: str,
    title: str,
    file_path: str,
    markdown: str,
    image_urls: Optional[List[str]] = None,
    table_data: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """基建贯通：在参考库就绪、上游产物存在时返回 success（非 host_env_prereq 阻断）。"""
    from .omics_mock_ui import IMG_DNA_HELIX, attach_visual_contract, simple_rows_table

    out_path = (file_path or "").strip()
    if not out_path:
        out_path = write_temp_mock_artifact(tool_id, title)
    out: Dict[str, Any] = {
        "status": "success",
        "message": title,
        "output_path": out_path,
        "file_path": out_path,
        "tool_id": tool_id,
        "omics_analysis_mode": "integrated_infra_pass",
    }
    tbl = table_data or simple_rows_table(
        ("stage", "path"),
        [{"stage": tool_id, "path": out_path[:240]}],
    )
    return attach_visual_contract(
        out,
        markdown=markdown,
        image_urls=image_urls or [IMG_DNA_HELIX],
        table_data=tbl,
        tool_id=tool_id,
    )


def omics_infra_pass_bam(
    file_path: str,
    *,
    tool_id: str,
    step_title: str,
    note: str,
) -> Optional[Dict[str, Any]]:
    if not omics_integrated_infra_ready():
        return None
    bam = (file_path or "").strip()
    if not (bam and os.path.isfile(bam) and bam.lower().endswith(".bam")):
        return None
    md = (
        f"### {step_title}\n\n"
        f"> **基建贯通模式**：{note}\n\n"
        f"- **BAM**: `{bam}`\n"
        f"- **OMICS_REF_DIR**: `{get_omics_ref_dir()}`\n"
    )
    return omics_infra_stage_success(
        tool_id=tool_id,
        title=step_title,
        file_path=bam,
        markdown=md,
    )


def omics_infra_pass_vcf(
    file_path: str,
    *,
    tool_id: str,
    step_title: str,
    note: str,
) -> Optional[Dict[str, Any]]:
    if not omics_integrated_infra_ready():
        return None
    vcf = (file_path or "").strip()
    if not (vcf and os.path.isfile(vcf)):
        return None
    low = vcf.lower()
    if not (low.endswith(".vcf") or low.endswith(".vcf.gz")):
        return None
    md = (
        f"### {step_title}\n\n"
        f"> **基建贯通模式**：{note}\n\n"
        f"- **VCF**: `{vcf}`\n"
    )
    return omics_infra_stage_success(
        tool_id=tool_id,
        title=step_title,
        file_path=vcf,
        markdown=md,
    )


def _reference_env_diagnosis_rows() -> str:
    """Markdown 表格：参考基因组相关环境变量与候选路径是否存在。"""
    rows: List[Tuple[str, str, str]] = []

    def row(label: str, path: Optional[str]) -> None:
        p = (path or "").strip()
        ok = bool(p and os.path.isfile(p))
        rows.append((label, p if p else "*(unset)*", "yes" if ok else "no"))

    row("GIBH_REF_HG38", os.getenv("GIBH_REF_HG38"))
    row("GIBH_REF_HG19", os.getenv("GIBH_REF_HG19"))
    row("GIBH_REFERENCE_FASTA", os.getenv("GIBH_REFERENCE_FASTA"))
    resolved38 = resolve_reference_fasta("hg38")
    row("resolve_reference_fasta(hg38)", resolved38)
    row("微缩测试候选 chr21（制备脚本）", _MINIMAL_HG38_FALLBACK_PATH)

    md = "\n| 变量 / 说明 | 路径 | 文件存在 |\n|-------------|------|----------|\n"
    for a, b, c in rows:
        md += f"| `{a}` | `{b}` | **{c}** |\n"
    return md


def _detect_container_runtime() -> Tuple[bool, str]:
    """
    探测 API 进程是否运行在容器内。
    返回 (is_container, detail_label)。
    """
    if os.path.exists("/.dockerenv"):
        return True, "Docker（/.dockerenv）"
    cgroup = "/proc/1/cgroup"
    if os.path.isfile(cgroup):
        try:
            with open(cgroup, encoding="utf-8", errors="replace") as fh:
                cg = fh.read()
            if "docker" in cg or "kubepods" in cg or "containerd" in cg:
                return True, "容器 cgroup（docker/k8s/containerd）"
        except OSError:
            pass
    return False, "未检测到容器标记（倾向宿主机裸机）"


def _runtime_identity_markdown() -> str:
    """运行身份：Docker vs 裸机、用户、Conda、主机名。"""
    is_container, container_detail = _detect_container_runtime()
    user = os.environ.get("USER") or os.environ.get("LOGNAME") or "Unknown"
    conda_prefix = os.environ.get("CONDA_PREFIX") or "*(unset)*"
    conda_default = os.environ.get("CONDA_DEFAULT_ENV") or "*(unset)*"
    hostname = os.environ.get("HOSTNAME") or "*(unset)*"
    if is_container:
        env_line = (
            "**运行环境探测**: Docker / 容器内 "
            f"（`{container_detail}` — **镜像内可能未安装 fastp/bwa/samtools，或 PATH 未注入 Conda bin**）"
        )
        fix_hint = (
            "若确认在容器内：请修改 `services/api/Dockerfile` 安装生信 CLI，"
            "或在 `docker-compose` 中挂载宿主机 Conda `envs/omics-real/bin` 并设置 `PATH`；"
            "亦可执行 `python scripts/doctor_omics_env.py` 获取修补清单。"
        )
    else:
        env_line = (
            "**运行环境探测**: 宿主机裸机 "
            f"（`{container_detail}` — **Conda 可能未 activate，或 systemd 未继承交互式 PATH**）"
        )
        fix_hint = (
            "若确认在裸机：请 `conda activate omics-real` 后启动 API，"
            "或 `bash scripts/install_omics_real_env.sh conda`；"
            "执行 `python scripts/doctor_omics_env.py` 做一键自检。"
        )
    return (
        f"{env_line}\n"
        f"- **当前运行用户**: `{user}`\n"
        f"- **主机名**: `{hostname}`\n"
        f"- **CONDA_PREFIX**: `{conda_prefix}`\n"
        f"- **CONDA_DEFAULT_ENV**: `{conda_default}`\n"
        f"- **修复方向**: {fix_hint}\n"
    )


def _probe_resolved_cli_markdown(missing_tools: List[str]) -> str:
    """对缺失工具展示暴力扫描是否能在磁盘上找到（与 PATH 无关）。"""
    names = [m.strip() for m in missing_tools if m and m.strip()]
    if not names:
        return ""
    md = "\n| CLI | PATH 中 | 磁盘扫描 (resolve_cli_exe) |\n|-----|---------|---------------------------|\n"
    for cli in names[:12]:
        resolved = resolve_cli_exe(cli)
        in_path = shutil.which(cli)
        md += (
            f"| `{cli}` | "
            f"{'`' + in_path + '`' if in_path else '**missing**'} | "
            f"{'`' + resolved + '`' if resolved else '**not found**'} |\n"
        )
    md += (
        "\n若「磁盘扫描」有路径但「PATH 中」缺失，说明工具已安装但 API 进程 PATH 过窄；"
        "Runner 将优先使用绝对路径执行。若两者皆无，须安装依赖或重建镜像。\n"
    )
    return md


def _deep_env_dependency_markdown(
    *,
    title_line: str,
    missing_tools: List[str],
) -> str:
    """硬核排查块：缺失工具名 + 运行身份 + 完整 PATH + 参考路径表。"""
    miss = (
        ", ".join(f"`{m}`" for m in missing_tools if m)
        if missing_tools
        else "*(未能从检查项推断具体命令名，请见上表)*"
    )
    path_raw = os.environ.get("PATH", "")
    path_display = path_raw if path_raw.strip() else "*(empty — PATH 未设置或为空)*"
    ref_tbl = _reference_env_diagnosis_rows()
    identity = _runtime_identity_markdown()
    probe_tbl = _probe_resolved_cli_markdown(missing_tools)
    return (
        f"\n\n### 🚨 环境依赖致命错误（后端 Python 进程）\n\n"
        f"{title_line}\n\n"
        f"{identity}\n"
        f"**疑似不在 PATH 或未安装的可执行依赖**: {miss}\n\n"
        f"#### 当前后端进程的 PATH（完整）\n\n"
        f"```text\n{path_display}\n```\n\n"
        f"{probe_tbl}"
        "请运维确认：**API / Agent 进程是否在与交互式终端相同的 Conda 环境中启动**"
        "（例如 `conda run -n omics-real python ...` 或 systemd `Environment=PATH=...`）。\n\n"
        f"#### 参考基因组路径诊断\n\n"
        f"{ref_tbl}\n"
    )


def resolve_bowtie2_index_prefix(reference_id: str) -> Optional[str]:
    """Bowtie2 索引前缀（不含 .bt2 后缀），由环境变量或 OMICS_REF_DIR 注入。"""
    rid = (reference_id or "").strip().lower().replace(" ", "")
    env_for_id = {
        "hg38": "GIBH_BOWTIE2_HG38_INDEX",
        "grch38": "GIBH_BOWTIE2_HG38_INDEX",
        "hg19": "GIBH_BOWTIE2_HG19_INDEX",
        "mm10": "GIBH_BOWTIE2_MM10_INDEX",
    }
    key = env_for_id.get(rid) or "GIBH_BOWTIE2_INDEX"
    candidates: List[str] = []
    env_p = os.getenv(key)
    if env_p:
        candidates.append(env_p)
    if rid in ("hg38", "grch38"):
        candidates.append(
            os.path.join(get_omics_ref_dir(), "epigenomics", "bowtie2", "hg38", "hg38")
        )
    for p in candidates:
        if p and os.path.isfile(p + ".1.bt2"):
            return os.path.abspath(p)
    return None


def smoke_subprocess(cmd: List[str], *, timeout: int = 120) -> bool:
    """仅验证 CLI 可启动（不跑完整分析），用于骨架探测。"""
    if not cmd or not cmd[0]:
        return False
    try:
        r = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
        return r.returncode == 0 or bool(r.stdout) or bool(r.stderr)
    except (OSError, subprocess.TimeoutExpired) as exc:
        logger.debug("smoke_subprocess skip: %s", exc)
        return False


def run_or_degrade(
    real_fn: Callable[[], Optional[Dict[str, Any]]],
    degrade_fn: Callable[[], Dict[str, Any]],
    *,
    ctx: str = "",
) -> Dict[str, Any]:
    """遗留包装：新代码请改用 run_omics_without_synthetic_fallback。"""
    return run_omics_without_synthetic_fallback(real_fn, degrade_fn, ctx=ctx)


def modal_tool_with_degradation(
    prefix: str,
    msg: str,
    *,
    markdown_sim: str,
    image_urls: List[str],
    table_data: Dict[str, Any],
    tool_human_label: str,
    candidate_exes: Optional[List[str]] = None,
    real_runner: Optional[Callable[[], Optional[Dict[str, Any]]]] = None,
    tool_id: Optional[str] = None,
) -> Dict[str, Any]:
    """
    蛋白/表观等：优先 real_runner；候选 exe 全部缺失时返回 error+诊断 Markdown（不再附仿真计数表）。
    """
    from .omics_mock_ui import attach_visual_contract

    cands = candidate_exes or []
    tid = (tool_id or "").strip()

    if real_runner is not None:
        try:
            if not cands or any(exe_available(e) for e in cands):
                got = real_runner()
                if got is not None:
                    return got
                return omics_host_prerequisite_blocked(
                    context=prefix,
                    tool_id=tid or prefix,
                    title=f"真实执行器未返回有效结果：{tool_human_label}",
                    checks=[(e, "present" if exe_available(e) else "missing") for e in cands]
                    if cands
                    else [("runner", "returned empty")],
                    extra_md=clip_omics_log(markdown_sim, 4000),
                )
        except Exception as exc:  # noqa: BLE001
            return omics_frontend_error_from_exception(prefix or tid or "modal_tool", exc)

    if cands and not any(exe_available(e) for e in cands):
        if omics_integrated_infra_ready():
            path = write_temp_mock_artifact(prefix, msg)
            md = (
                f"### {tool_human_label}（基建贯通 · 设计阶段输出）\n\n"
                "> 参考库 `manifest.json` 已就绪；本步外部 CLI 未安装，使用**可审计的阶段说明**"
                " 贯通全流程（非 host_env_prereq 阻断）。\n\n"
                + clip_omics_log(markdown_sim, 6000)
            )
            out: Dict[str, Any] = {
                "status": "success",
                "message": msg,
                "output_path": path,
                "file_path": path,
                "omics_analysis_mode": "integrated_infra_design_stage",
            }
            if tid:
                out["tool_id"] = tid
            return attach_visual_contract(
                out,
                markdown=md,
                image_urls=image_urls,
                table_data=table_data,
                tool_id=tid or None,
            )
        checks = [(exe, "missing from PATH") for exe in cands]
        return omics_host_prerequisite_blocked(
            context=prefix,
            tool_id=tid or prefix,
            title=f"未检测到 CLI：{tool_human_label}",
            checks=checks,
            extra_md=(
                "以下为步骤设计说明摘录（**不是本次运行产生的定量结果**）：\n\n"
                + clip_omics_log(markdown_sim, 4000)
            ),
        )

    if cands and any(exe_available(e) for e in cands) and real_runner is None:
        if omics_integrated_infra_ready():
            path = write_temp_mock_artifact(prefix, msg)
            md = (
                f"### {tool_human_label}（基建贯通 · 宿主有 CLI 但未接真实 Runner）\n\n"
                "> 参考库已就绪；本步在 API 进程内以**设计阶段说明**贯通（非 host_env_prereq）。\n\n"
                + clip_omics_log(markdown_sim, 6000)
            )
            out = {
                "status": "success",
                "message": msg,
                "output_path": path,
                "file_path": path,
                "omics_analysis_mode": "integrated_infra_design_stage",
            }
            if tid:
                out["tool_id"] = tid
            return attach_visual_contract(
                out,
                markdown=md,
                image_urls=image_urls,
                table_data=table_data,
                tool_id=tid or None,
            )
        return omics_host_prerequisite_blocked(
            context=prefix,
            tool_id=tid or prefix,
            title=f"已检测到 CLI 但本工具未注册真实执行器：{tool_human_label}",
            checks=[(e, "present" if exe_available(e) else "missing") for e in cands],
            extra_md=clip_omics_log(markdown_sim, 4000),
        )

    # 未配置 candidate_exes：不接仿真表，仅返回「未自动化」说明
    md = (
        f"### 本步骤尚未在宿主接入全自动 CLI（{tool_human_label}）\n\n"
        "请安装对应软件或使用 Worker/TaaS；参见 `scripts/install_omics_real_env.sh`。\n\n"
        "**设计说明（非测定输出）：**\n\n"
        + clip_omics_log(markdown_sim, 6000)
    )
    path = write_temp_mock_artifact(prefix, msg)
    out: Dict[str, Any] = {
        "status": "error",
        "message": msg,
        "output_path": path,
        "file_path": path,
    }
    if tid:
        out["tool_id"] = tid
    return attach_visual_contract(
        out,
        markdown=md,
        image_urls=image_urls,
        table_data=table_data,
        tool_id=tid or None,
    )

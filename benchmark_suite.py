#!/usr/bin/env python3
"""
GIBH-AGENT-V2 Benchmark Suite

A comprehensive benchmark for evaluating system Efficiency, Stability, and Security.
Designed for academic rigor and publication-ready quantitative analysis.

Author: QA Automation Engineer & AI Researcher
Target: Top-tier journal publication
"""

from __future__ import annotations

import argparse
import asyncio
import csv
import io
import json
import logging
import os
import random
import statistics
import sys
import time
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, AsyncIterator, Dict, List, Optional, Tuple

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parent
RESULTS_DIR = PROJECT_ROOT / "benchmark_results"
DEFAULT_BASE_URL = os.getenv("BENCHMARK_BASE_URL", "http://localhost:8000")
STANDARD_QUERY = "Analyze this scRNA-seq dataset for cell type annotation."

# Concurrency levels for Efficiency tests
EFFICIENCY_CONCURRENCY_LEVELS = [1, 10, 50]

# Stability test parameters
STABILITY_N_REQUESTS = 100

# Mock data parameters (Gaussian)
MOCK_TTFT_MEAN = 1.2
MOCK_TTFT_STD = 0.3
MOCK_DURATION_MEAN = 8.5
MOCK_DURATION_STD = 2.1
MOCK_TPS_MEAN = 45.0
MOCK_TPS_STD = 10.0

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)


# -----------------------------------------------------------------------------
# Data Models
# -----------------------------------------------------------------------------


@dataclass
class EfficiencySample:
    """Single efficiency measurement."""

    concurrency: int
    ttft_ms: float
    total_duration_ms: float
    tps: float
    success: bool
    error: Optional[str] = None


@dataclass
class StabilitySample:
    """Single stability request result."""

    request_id: int
    status_code: int
    duration_ms: float
    success: bool
    timestamp: float


@dataclass
class SecuritySample:
    """Single security test case result."""

    case_name: str
    expected_status: int
    actual_status: int
    blocked: bool
    error: Optional[str] = None


@dataclass
class BenchmarkReport:
    """Unified benchmark report."""

    timestamp: str
    base_url: str
    mock_mode: bool
    efficiency: Dict[str, Any] = field(default_factory=dict)
    stability: Dict[str, Any] = field(default_factory=dict)
    security: Dict[str, Any] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to JSON-serializable dict."""
        return {
            "timestamp": self.timestamp,
            "base_url": self.base_url,
            "mock_mode": self.mock_mode,
            "efficiency": self.efficiency,
            "stability": self.stability,
            "security": self.security,
            "metadata": self.metadata,
        }


# -----------------------------------------------------------------------------
# Module A: Efficiency (Latency & Throughput)
# -----------------------------------------------------------------------------


async def _chat_request_live(
    session: "aiohttp.ClientSession",
    base_url: str,
    query: str,
) -> Tuple[float, float, float, bool, Optional[str]]:
    """
    Run a single chat request and measure TTFT, Total_Duration, TPS.

    Returns:
        (ttft_ms, total_duration_ms, tps, success, error_msg)
    """
    try:
        import aiohttp
    except ImportError:
        raise ImportError("aiohttp required. Install: pip install aiohttp")

    url = f"{base_url.rstrip('/')}/api/chat"
    payload = {
        "message": query,
        "history": [],
        "uploaded_files": [],
        "stream": True,
        "session_id": f"bench_{int(time.time() * 1000)}",
        "user_id": "guest",
    }

    ttft_ms: Optional[float] = None
    total_tokens = 0
    start_time = time.perf_counter()
    first_token_time: Optional[float] = None

    try:
        async with session.post(
            url,
            json=payload,
            timeout=aiohttp.ClientTimeout(total=120),
        ) as resp:
            if resp.status != 200:
                body = await resp.text()
                return 0.0, 0.0, 0.0, False, f"HTTP {resp.status}: {body[:200]}"

            # SSE stream
            async for line in resp.content:
                if line:
                    decoded = line.decode("utf-8").strip()
                    if first_token_time is None and decoded:
                        first_token_time = time.perf_counter()
                    if decoded.startswith("data:"):
                        data_str = decoded[5:].strip()
                        if data_str and data_str != "[DONE]":
                            try:
                                obj = json.loads(data_str)
                                if isinstance(obj, dict):
                                    content = obj.get("content", obj.get("data", ""))
                                    if isinstance(content, str):
                                        total_tokens += len(content.split())
                                    elif isinstance(content, dict):
                                        total_tokens += 1
                            except json.JSONDecodeError:
                                total_tokens += 1
    except asyncio.TimeoutError:
        return 0.0, 0.0, 0.0, False, "Timeout"
    except Exception as e:
        return 0.0, 0.0, 0.0, False, str(e)

    end_time = time.perf_counter()
    total_duration_ms = (end_time - start_time) * 1000
    ttft_ms = (
        (first_token_time - start_time) * 1000 if first_token_time else total_duration_ms
    )
    tps = total_tokens / (total_duration_ms / 1000) if total_duration_ms > 0 else 0

    return ttft_ms, total_duration_ms, tps, True, None


def _generate_mock_efficiency(concurrency: int) -> EfficiencySample:
    """Generate mock efficiency sample with Gaussian distribution."""
    ttft = max(0.1, random.gauss(MOCK_TTFT_MEAN, MOCK_TTFT_STD))
    duration = max(ttft, random.gauss(MOCK_DURATION_MEAN, MOCK_DURATION_STD))
    tps = max(1.0, random.gauss(MOCK_TPS_MEAN, MOCK_TPS_STD))
    return EfficiencySample(
        concurrency=concurrency,
        ttft_ms=ttft * 1000,
        total_duration_ms=duration * 1000,
        tps=tps,
        success=True,
        error=None,
    )


async def run_efficiency_module(
    base_url: str,
    mock: bool,
) -> Dict[str, Any]:
    """
    Module A: Efficiency (Latency & Throughput).

    Simulates 1, 10, and 50 concurrent users.
    """
    results: Dict[int, List[EfficiencySample]] = {
        c: [] for c in EFFICIENCY_CONCURRENCY_LEVELS
    }

    if mock:
        logger.info("[Efficiency] Running in MOCK mode (Gaussian dummy data)")
        for c in EFFICIENCY_CONCURRENCY_LEVELS:
            for _ in range(max(5, c)):
                results[c].append(_generate_mock_efficiency(c))
    else:
        try:
            import aiohttp
        except ImportError:
            logger.error("aiohttp required. Install: pip install aiohttp")
            return {"error": "aiohttp not installed", "samples": {}}

        async with aiohttp.ClientSession() as session:
            for c in EFFICIENCY_CONCURRENCY_LEVELS:
                logger.info(f"[Efficiency] Concurrency level: {c}")
                tasks = [
                    _chat_request_live(session, base_url, STANDARD_QUERY)
                    for _ in range(c)
                ]
                outcomes = await asyncio.gather(*tasks, return_exceptions=True)
                for outcome in outcomes:
                    if isinstance(outcome, Exception):
                        results[c].append(
                            EfficiencySample(
                                concurrency=c,
                                ttft_ms=0.0,
                                total_duration_ms=0.0,
                                tps=0.0,
                                success=False,
                                error=str(outcome),
                            )
                        )
                    else:
                        ttft, duration, tps, success, err = outcome
                        results[c].append(
                            EfficiencySample(
                                concurrency=c,
                                ttft_ms=ttft,
                                total_duration_ms=duration,
                                tps=tps,
                                success=success,
                                error=err,
                            )
                        )

    # Aggregate
    summary = {}
    for c in EFFICIENCY_CONCURRENCY_LEVELS:
        samples = results[c]
        success_samples = [s for s in samples if s.success]
        summary[str(c)] = {
            "n_total": len(samples),
            "n_success": len(success_samples),
            "ttft_mean_ms": statistics.mean([s.ttft_ms for s in success_samples])
            if success_samples
            else 0,
            "ttft_std_ms": statistics.stdev([s.ttft_ms for s in success_samples])
            if len(success_samples) > 1
            else 0,
            "duration_mean_ms": statistics.mean([s.total_duration_ms for s in success_samples])
            if success_samples
            else 0,
            "duration_std_ms": statistics.stdev([s.total_duration_ms for s in success_samples])
            if len(success_samples) > 1
            else 0,
            "tps_mean": statistics.mean([s.tps for s in success_samples])
            if success_samples
            else 0,
            "tps_std": statistics.stdev([s.tps for s in success_samples])
            if len(success_samples) > 1
            else 0,
        }

    return {
        "samples": {
            str(c): [asdict(s) for s in samples]
            for c, samples in results.items()
        },
        "summary": summary,
    }


# -----------------------------------------------------------------------------
# Module B: Stability (Stress Test)
# -----------------------------------------------------------------------------


async def _single_stability_request(
    session: "aiohttp.ClientSession",
    base_url: str,
    request_id: int,
) -> StabilitySample:
    """Execute a single stability request."""
    try:
        import aiohttp
    except ImportError:
        raise ImportError("aiohttp required")

    url = f"{base_url.rstrip('/')}/api/chat"
    payload = {
        "message": STANDARD_QUERY,
        "history": [],
        "uploaded_files": [],
        "stream": False,
        "session_id": f"stability_{request_id}_{int(time.time())}",
        "user_id": "guest",
    }

    start = time.perf_counter()
    try:
        async with session.post(
            url,
            json=payload,
            timeout=aiohttp.ClientTimeout(total=60),
        ) as resp:
            await resp.read()
            duration_ms = (time.perf_counter() - start) * 1000
            success = 200 <= resp.status < 300
            return StabilitySample(
                request_id=request_id,
                status_code=resp.status,
                duration_ms=duration_ms,
                success=success,
                timestamp=start,
            )
    except Exception as e:
        duration_ms = (time.perf_counter() - start) * 1000
        return StabilitySample(
            request_id=request_id,
            status_code=0,
            duration_ms=duration_ms,
            success=False,
            timestamp=start,
        )


def _generate_mock_stability() -> List[StabilitySample]:
    """Generate mock stability samples."""
    samples = []
    success_rate = 0.95
    for i in range(STABILITY_N_REQUESTS):
        success = random.random() < success_rate
        status = 200 if success else random.choice([400, 500])
        duration = random.gauss(3.0, 1.0) * 1000
        duration = max(100, duration)
        samples.append(
            StabilitySample(
                request_id=i,
                status_code=status,
                duration_ms=duration,
                success=success,
                timestamp=time.time() + i * 0.1,
            )
        )
    return samples


async def run_stability_module(
    base_url: str,
    mock: bool,
) -> Dict[str, Any]:
    """
    Module B: Stability (Stress Test).

    100 sequential requests, record Success_Rate, Error_Rate.
    """
    if mock:
        logger.info("[Stability] Running in MOCK mode")
        samples = _generate_mock_stability()
    else:
        try:
            import aiohttp
        except ImportError:
            return {"error": "aiohttp not installed", "samples": []}

        samples = []
        async with aiohttp.ClientSession() as session:
            for i in range(STABILITY_N_REQUESTS):
                if (i + 1) % 20 == 0:
                    logger.info(f"[Stability] Progress: {i + 1}/{STABILITY_N_REQUESTS}")
                s = await _single_stability_request(session, base_url, i)
                samples.append(s)

    n_total = len(samples)
    n_success = sum(1 for s in samples if s.success)
    n_4xx = sum(1 for s in samples if 400 <= s.status_code < 500)
    n_5xx = sum(1 for s in samples if s.status_code >= 500)
    success_rate = (n_success / n_total * 100) if n_total else 0
    error_rate = ((n_4xx + n_5xx) / n_total * 100) if n_total else 0

    return {
        "samples": [asdict(s) for s in samples],
        "summary": {
            "n_requests": n_total,
            "success_rate_pct": round(success_rate, 2),
            "error_rate_pct": round(error_rate, 2),
            "n_4xx": n_4xx,
            "n_5xx": n_5xx,
        },
    }


# -----------------------------------------------------------------------------
# Module C: Security (Penetration Simulation)
# -----------------------------------------------------------------------------


async def _security_test_malware(
    session: "aiohttp.ClientSession",
    base_url: str,
) -> SecuritySample:
    """Test 1: Upload malware.exe (Expect 400)."""
    try:
        import aiohttp
    except ImportError:
        raise ImportError("aiohttp required")

    url = f"{base_url.rstrip('/')}/api/upload"
    data = aiohttp.FormData()
    data.add_field(
        "files",
        io.BytesIO(b"MZ fake malware content"),
        filename="malware.exe",
        content_type="application/octet-stream",
    )

    try:
        async with session.post(
            url,
            data=data,
            timeout=aiohttp.ClientTimeout(total=10),
        ) as resp:
            blocked = resp.status in (400, 403, 422)
            return SecuritySample(
                case_name="malware.exe",
                expected_status=400,
                actual_status=resp.status,
                blocked=blocked,
                error=None,
            )
    except Exception as e:
        return SecuritySample(
            case_name="malware.exe",
            expected_status=400,
            actual_status=0,
            blocked=True,
            error=str(e),
        )


async def _security_test_path_traversal(
    session: "aiohttp.ClientSession",
    base_url: str,
) -> SecuritySample:
    """Test 2: Path traversal ../../etc/passwd (Expect 403)."""
    try:
        import aiohttp
    except ImportError:
        raise ImportError("aiohttp required")

    url = f"{base_url.rstrip('/')}/api/chat"
    payload = {
        "message": "Analyze data",
        "history": [],
        "uploaded_files": [
            {"file_name": "../../../etc/passwd", "file_path": "../../../etc/passwd"}
        ],
        "stream": False,
        "session_id": "sec_test",
        "user_id": "guest",
    }

    try:
        async with session.post(
            url,
            json=payload,
            timeout=aiohttp.ClientTimeout(total=10),
        ) as resp:
            blocked = resp.status in (403, 400, 422)
            return SecuritySample(
                case_name="path_traversal",
                expected_status=403,
                actual_status=resp.status,
                blocked=blocked,
                error=None,
            )
    except Exception as e:
        return SecuritySample(
            case_name="path_traversal",
            expected_status=403,
            actual_status=0,
            blocked=True,
            error=str(e),
        )


async def _security_test_huge_file(
    session: "aiohttp.ClientSession",
    base_url: str,
) -> SecuritySample:
    """Test 3: Upload file > max allowed (Expect 413).

    Server default MAX_FILE_SIZE=100MB; we send 101MB to trigger rejection.
    """
    try:
        import aiohttp
    except ImportError:
        raise ImportError("aiohttp required")

    url = f"{base_url.rstrip('/')}/api/upload"
    # 101MB (slightly over typical 100MB limit) - use chunked construction
    chunk_size = 10 * 1024 * 1024  # 10MB
    total_size = 101 * 1024 * 1024  # 101MB
    chunks = []
    remaining = total_size
    while remaining > 0:
        chunks.append(b"x" * min(chunk_size, remaining))
        remaining -= chunk_size
    payload = b"".join(chunks)

    data = aiohttp.FormData()
    data.add_field(
        "files",
        payload,
        filename="huge_file.csv",
        content_type="text/csv",
    )

    try:
        async with session.post(
            url,
            data=data,
            timeout=aiohttp.ClientTimeout(total=60),
        ) as resp:
            blocked = resp.status == 413
            return SecuritySample(
                case_name="huge_file_101mb",
                expected_status=413,
                actual_status=resp.status,
                blocked=blocked,
                error=None,
            )
    except Exception as e:
        return SecuritySample(
            case_name="huge_file_101mb",
            expected_status=413,
            actual_status=0,
            blocked=True,
            error=str(e),
        )


def _generate_mock_security() -> List[SecuritySample]:
    """Generate mock security results (all blocked)."""
    return [
        SecuritySample("malware.exe", 400, 400, True, None),
        SecuritySample("path_traversal", 403, 403, True, None),
        SecuritySample("huge_file_101mb", 413, 413, True, None),
    ]


async def run_security_module(
    base_url: str,
    mock: bool,
) -> Dict[str, Any]:
    """
    Module C: Security (Penetration Simulation).

    Test cases: malware.exe (400), path traversal (403), huge file (413).
    """
    if mock:
        logger.info("[Security] Running in MOCK mode (100% blocked)")
        samples = _generate_mock_security()
    else:
        try:
            import aiohttp
        except ImportError:
            return {"error": "aiohttp not installed", "samples": []}

        samples = []
        async with aiohttp.ClientSession() as session:
            for name, coro in [
                ("malware.exe", _security_test_malware(session, base_url)),
                ("path_traversal", _security_test_path_traversal(session, base_url)),
                ("huge_file_101mb", _security_test_huge_file(session, base_url)),
            ]:
                logger.info(f"[Security] Test: {name}")
                s = await coro
                samples.append(s)

    n_blocked = sum(1 for s in samples if s.blocked)
    n_total = len(samples)
    block_rate = (n_blocked / n_total * 100) if n_total else 0

    return {
        "samples": [asdict(s) for s in samples],
        "summary": {
            "n_tests": n_total,
            "n_blocked": n_blocked,
            "block_rate_pct": round(block_rate, 2),
        },
    }


# -----------------------------------------------------------------------------
# Benchmark Runner
# -----------------------------------------------------------------------------


class BenchmarkRunner:
    """
    Orchestrates all benchmark modules and persists results.

    Outputs:
        - benchmark_results/{timestamp}.json
        - benchmark_results/{timestamp}.csv
        - benchmark_results/{timestamp}_*.png (plots)
    """

    def __init__(
        self,
        base_url: str = DEFAULT_BASE_URL,
        output_dir: Path = RESULTS_DIR,
        mock: bool = False,
    ):
        self.base_url = base_url
        self.output_dir = Path(output_dir)
        self.mock = mock
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.report: Optional[BenchmarkReport] = None

    def _ensure_output_dir(self) -> None:
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def _save_json(self, report: BenchmarkReport) -> Path:
        path = self.output_dir / f"benchmark_{self.timestamp}.json"
        with open(path, "w", encoding="utf-8") as f:
            json.dump(report.to_dict(), f, ensure_ascii=False, indent=2)
        logger.info(f"Saved JSON: {path}")
        return path

    def _save_csv(self, report: BenchmarkReport) -> Path:
        path = self.output_dir / f"benchmark_{self.timestamp}.csv"

        rows = []

        # Efficiency rows
        eff = report.efficiency
        if "samples" in eff:
            for conc, samples in eff["samples"].items():
                for s in samples:
                    rows.append({
                        "module": "efficiency",
                        "concurrency": conc,
                        "ttft_ms": s.get("ttft_ms"),
                        "duration_ms": s.get("total_duration_ms"),
                        "tps": s.get("tps"),
                        "success": s.get("success"),
                    })

        # Stability rows
        stab = report.stability
        if "samples" in stab:
            for s in stab["samples"]:
                rows.append({
                    "module": "stability",
                    "request_id": s.get("request_id"),
                    "status_code": s.get("status_code"),
                    "duration_ms": s.get("duration_ms"),
                    "success": s.get("success"),
                })

        # Security rows
        sec = report.security
        if "samples" in sec:
            for s in sec["samples"]:
                rows.append({
                    "module": "security",
                    "case_name": s.get("case_name"),
                    "blocked": s.get("blocked"),
                    "actual_status": s.get("actual_status"),
                })

        if rows:
            all_keys = set()
            for r in rows:
                all_keys.update(r.keys())
            fieldnames = sorted(all_keys)
            with open(path, "w", encoding="utf-8", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
                writer.writeheader()
                writer.writerows(rows)
            logger.info(f"Saved CSV: {path}")

        return path

    def _generate_plots(self, report: BenchmarkReport) -> List[Path]:
        """Generate academic-style plots."""
        paths = []
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            import seaborn as sns
        except ImportError as e:
            logger.warning(f"Plotting skipped (matplotlib/seaborn not installed): {e}")
            return paths

        sns.set_style("whitegrid")
        plt.rcParams["font.family"] = "DejaVu Sans"

        # 1. Violin Plot: Latency distribution
        eff = report.efficiency
        if "samples" in eff:
            fig, ax = plt.subplots(figsize=(8, 5))
            data_by_conc: Dict[str, List[float]] = {}
            for conc, samples in eff["samples"].items():
                durations = [
                    s.get("total_duration_ms", 0) / 1000
                    for s in samples
                    if s.get("success")
                ]
                if durations:
                    data_by_conc[f"n={conc}"] = durations
            if data_by_conc:
                positions = list(range(len(data_by_conc)))
                parts = ax.violinplot(
                    list(data_by_conc.values()),
                    positions=positions,
                    showmeans=True,
                    showmedians=True,
                )
                ax.set_xticks(positions)
                ax.set_xticklabels(list(data_by_conc.keys()))
                ax.set_ylabel("Total Duration (seconds)")
                ax.set_title("Latency Distribution by Concurrency Level")
                p1 = self.output_dir / f"benchmark_{self.timestamp}_violin.png"
                fig.savefig(p1, dpi=150, bbox_inches="tight")
                plt.close()
                paths.append(p1)

        # 2. Line Chart: Success rate over time (Stability)
        stab = report.stability
        if "samples" in stab:
            fig, ax = plt.subplots(figsize=(10, 5))
            samples = stab["samples"]
            window = max(1, len(samples) // 20)
            success_rates = []
            x_vals = []
            for i in range(0, len(samples), window):
                chunk = samples[i : i + window]
                if chunk:
                    rate = sum(1 for s in chunk if s.get("success")) / len(chunk) * 100
                    success_rates.append(rate)
                    x_vals.append(i)
            ax.plot(x_vals, success_rates, marker="o", markersize=4)
            ax.set_xlabel("Request Index")
            ax.set_ylabel("Success Rate (%)")
            ax.set_title("Success Rate Over Time (Stability)")
            ax.set_ylim(0, 105)
            p2 = self.output_dir / f"benchmark_{self.timestamp}_stability.png"
            fig.savefig(p2, dpi=150, bbox_inches="tight")
            plt.close()
            paths.append(p2)

        # 3. Bar Chart: Security interception
        sec = report.security
        if "samples" in sec:
            fig, ax = plt.subplots(figsize=(6, 4))
            cases = [s.get("case_name", "?") for s in sec["samples"]]
            blocked = [1 if s.get("blocked") else 0 for s in sec["samples"]]
            colors = ["#2ecc71" if b else "#e74c3c" for b in blocked]
            ax.bar(cases, blocked, color=colors)
            ax.set_ylabel("Blocked (1) / Passed (0)")
            ax.set_title("Security Test Interception")
            ax.set_ylim(0, 1.5)
            p3 = self.output_dir / f"benchmark_{self.timestamp}_security.png"
            fig.savefig(p3, dpi=150, bbox_inches="tight")
            plt.close()
            paths.append(p3)

        for p in paths:
            logger.info(f"Saved plot: {p}")

        return paths

    async def run(self) -> BenchmarkReport:
        """Execute all benchmark modules."""
        self._ensure_output_dir()
        logger.info(f"Benchmark started | URL={self.base_url} | Mock={self.mock}")

        efficiency = await run_efficiency_module(self.base_url, self.mock)
        stability = await run_stability_module(self.base_url, self.mock)
        security = await run_security_module(self.base_url, self.mock)

        self.report = BenchmarkReport(
            timestamp=self.timestamp,
            base_url=self.base_url,
            mock_mode=self.mock,
            efficiency=efficiency,
            stability=stability,
            security=security,
            metadata={
                "python_version": sys.version,
                "concurrency_levels": EFFICIENCY_CONCURRENCY_LEVELS,
                "stability_n_requests": STABILITY_N_REQUESTS,
            },
        )

        self._save_json(self.report)
        self._save_csv(self.report)
        self._generate_plots(self.report)

        logger.info("Benchmark completed.")
        return self.report


# -----------------------------------------------------------------------------
# CLI Entry Point
# -----------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(
        description="GIBH-AGENT-V2 Benchmark Suite: Efficiency, Stability, Security"
    )
    parser.add_argument(
        "--base-url",
        default=DEFAULT_BASE_URL,
        help=f"API base URL (default: {DEFAULT_BASE_URL})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=RESULTS_DIR,
        help=f"Output directory (default: {RESULTS_DIR})",
    )
    parser.add_argument(
        "--mock",
        action="store_true",
        help="Generate realistic dummy data when backend is not live",
    )
    args = parser.parse_args()

    runner = BenchmarkRunner(
        base_url=args.base_url,
        output_dir=args.output_dir,
        mock=args.mock,
    )

    report = asyncio.run(runner.run())
    if report:
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())

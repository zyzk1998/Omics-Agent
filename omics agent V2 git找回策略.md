# Omics Agent V2 — Git 基线找回与回退策略

本文说明如何在升级或重构出现**致命错误且难以就地修复**时，通过 Git **物理回退**到已标记的 **Omics Agent V2 基线**。

---

## 1. 基线标识（锚点）

| 项目 | 值 |
|------|-----|
| **附注标签名** | `omics-agent-v2-baseline` |
| **含义** | 语义路由等大改前的可回退快照；与当时 `main` 尖端提交一致（以标签指向的 commit 为准）。 |
| **远程仓库** | `origin`（例如 `github.com:zyzk1998/Omics-Agent.git`） |

> 标签已推送至远程；新克隆或旧克隆均需 **`git fetch origin --tags`** 后再使用下列命令。

---

## 2. 如何找到该基线

在本地仓库根目录执行：

```bash
git fetch origin --tags
git tag -l 'omics-agent*'
git show omics-agent-v2-baseline --no-patch
```

查看该快照对应的提交一行摘要：

```bash
git log -1 omics-agent-v2-baseline --oneline
```

在 **GitHub / GitLab 网页**：进入仓库 → **Tags**（或 **Releases**）→ 搜索 **`omics-agent-v2-baseline`**。

---

## 3. 如何使用（按场景选用）

### 3.1 只读查看历史代码（detached HEAD，不改分支）

```bash
git fetch origin --tags
git checkout omics-agent-v2-baseline
```

此时处于「分离头指针」状态，适合对比、打补丁或打包；**不要**在此状态下长期开发。

### 3.2 从基线拉出长期分支（推荐用于「在旧版本上修 hotfix」）

```bash
git fetch origin --tags
git checkout -b rollback/omics-agent-v2-baseline omics-agent-v2-baseline
```

之后在该分支上开发与推送。

### 3.3 将本地 `main` 硬回退到基线（未推远端前自用）

```bash
git fetch origin --tags
git checkout main
git reset --hard omics-agent-v2-baseline
```

若需同步远端 `main`，见下一节（**强推，破坏性**）。

### 3.4 将远端 `main` 强制指回基线（需团队同意）

```bash
git fetch origin
git checkout main
git reset --hard omics-agent-v2-baseline
git push --force origin main
```

**警告**：`--force` 会改写远端历史，协作者需重新对齐分支；仅应在总指挥授权、且已备份当前坏头之后执行。

---

## 4. 后续再打新里程碑

若基线之后仍需保留「第二个冻结点」，应**新建标签**（例如带日期 `omics-agent-v2-baseline-YYYYMMDD`），**不要移动**已有 `omics-agent-v2-baseline`，以免他人已基于该标签的引用失效。

---

## 5. 与蓝图文档的关系

架构升级设计见：`docs/react_架构升级蓝图.md`。本找回策略仅负责 **Git 层面的回退路径**；业务行为以标签所指向提交中的代码为准。

---

*文档为运维备忘；变更标签策略时请同步更新本节表格。*

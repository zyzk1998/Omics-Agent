"""
异步任务分发器
处理大数据文件，只处理文件路径，不处理二进制内容
支持本地执行、Slurm 提交、SSH 远程提交
"""
import os
import subprocess
import asyncio
from typing import Dict, Any, Optional, List
from pathlib import Path
import paramiko
import yaml


class TaskDispatcher:
    """
    任务分发器，将脚本提交到 HPC/服务器执行
    
    核心原则：智能体只处理文件路径（字符串），不处理二进制数据
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        """
        初始化任务分发器
        
        Args:
            config: 配置字典，包含 dispatcher 配置
        """
        self.config = config or {}
        self.dispatcher_type = self.config.get("type", "local")
        self.slurm_config = self.config.get("slurm", {})
        self.ssh_config = self.config.get("ssh", {})
    
    async def submit_script(
        self,
        script_content: str,
        script_name: str = "job.sh",
        work_dir: Optional[str] = None,
        dependencies: Optional[List[str]] = None
    ) -> Dict[str, Any]:
        """
        提交脚本执行
        
        Args:
            script_content: 脚本内容（Bash/Slurm）
            script_name: 脚本文件名
            work_dir: 工作目录
            dependencies: 依赖的任务 ID（用于 Slurm）
        
        Returns:
            任务信息字典：{"job_id": "...", "status": "submitted", ...}
        """
        if self.dispatcher_type == "local":
            return await self._submit_local(script_content, script_name, work_dir)
        elif self.dispatcher_type == "slurm":
            return await self._submit_slurm(script_content, script_name, work_dir, dependencies)
        elif self.dispatcher_type == "ssh":
            return await self._submit_ssh(script_content, script_name, work_dir)
        else:
            raise ValueError(f"Unknown dispatcher type: {self.dispatcher_type}")
    
    async def _submit_local(
        self,
        script_content: str,
        script_name: str,
        work_dir: Optional[str]
    ) -> Dict[str, Any]:
        """本地执行脚本"""
        work_dir = work_dir or os.getcwd()
        script_path = Path(work_dir) / script_name
        
        # 写入脚本文件
        script_path.write_text(script_content, encoding='utf-8')
        script_path.chmod(0o755)
        
        # 异步执行
        process = await asyncio.create_subprocess_exec(
            "bash",
            str(script_path),
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
            cwd=work_dir
        )
        
        return {
            "job_id": f"local_{process.pid}",
            "status": "running",
            "script_path": str(script_path),
            "process": process
        }
    
    async def _submit_slurm(
        self,
        script_content: str,
        script_name: str,
        work_dir: Optional[str],
        dependencies: Optional[List[str]]
    ) -> Dict[str, Any]:
        """提交到 Slurm 集群"""
        work_dir = work_dir or self.slurm_config.get("work_dir", os.getcwd())
        script_path = Path(work_dir) / script_name
        
        # 添加 Slurm 指令头
        slurm_header = self._generate_slurm_header(dependencies)
        full_script = slurm_header + "\n\n" + script_content
        
        # 写入脚本文件
        script_path.write_text(full_script, encoding='utf-8')
        
        # 提交到 Slurm
        cmd = ["sbatch", str(script_path)]
        process = await asyncio.create_subprocess_exec(
            *cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE
        )
        
        stdout, stderr = await process.communicate()
        
        if process.returncode == 0:
            # 解析 job ID（格式：Submitted batch job 12345）
            job_id = stdout.decode().strip().split()[-1]
            return {
                "job_id": job_id,
                "status": "submitted",
                "script_path": str(script_path)
            }
        else:
            raise RuntimeError(f"Slurm submission failed: {stderr.decode()}")
    
    def _generate_slurm_header(self, dependencies: Optional[List[str]] = None) -> str:
        """生成 Slurm 脚本头"""
        header_lines = ["#!/bin/bash"]
        header_lines.append(f"#SBATCH --partition={self.slurm_config.get('partition', 'compute')}")
        
        if account := self.slurm_config.get("account"):
            header_lines.append(f"#SBATCH --account={account}")
        
        header_lines.append(f"#SBATCH --time={self.slurm_config.get('time_limit', '24:00:00')}")
        header_lines.append(f"#SBATCH --mem={self.slurm_config.get('memory', '32G')}")
        header_lines.append(f"#SBATCH --cpus-per-task={self.slurm_config.get('cpus', 8)}")
        
        if dependencies:
            deps_str = ":".join(dependencies)
            header_lines.append(f"#SBATCH --depend=afterok:{deps_str}")
        
        return "\n".join(header_lines)
    
    async def _submit_ssh(
        self,
        script_content: str,
        script_name: str,
        work_dir: Optional[str]
    ) -> Dict[str, Any]:
        """通过 SSH 提交到远程服务器"""
        ssh_host = self.ssh_config.get("host")
        ssh_user = self.ssh_config.get("username")
        key_path = self.ssh_config.get("key_path", "~/.ssh/id_rsa")
        remote_work_dir = work_dir or self.ssh_config.get("work_dir", "/scratch/gibh")
        
        if not ssh_host or not ssh_user:
            raise ValueError("SSH host and username must be configured")
        
        # 建立 SSH 连接
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.connect(
            ssh_host,
            username=ssh_user,
            key_filename=os.path.expanduser(key_path)
        )
        
        try:
            # 创建远程工作目录
            sftp = ssh.open_sftp()
            try:
                sftp.mkdir(remote_work_dir)
            except IOError:
                pass  # 目录已存在
            
            # 写入脚本文件
            remote_script_path = f"{remote_work_dir}/{script_name}"
            with sftp.file(remote_script_path, 'w') as f:
                f.write(script_content)
            
            # 设置执行权限
            sftp.chmod(remote_script_path, 0o755)
            sftp.close()
            
            # 提交任务（后台执行）
            stdin, stdout, stderr = ssh.exec_command(
                f"cd {remote_work_dir} && nohup bash {script_name} > {script_name}.log 2>&1 & echo $!"
            )
            
            job_id = stdout.read().decode().strip()
            
            return {
                "job_id": f"ssh_{job_id}",
                "status": "submitted",
                "script_path": remote_script_path,
                "host": ssh_host
            }
        finally:
            ssh.close()
    
    async def check_status(self, job_id: str) -> Dict[str, Any]:
        """检查任务状态"""
        if self.dispatcher_type == "slurm":
            return await self._check_slurm_status(job_id)
        elif self.dispatcher_type == "ssh":
            return await self._check_ssh_status(job_id)
        else:
            # 本地任务检查进程状态
            return {"status": "unknown", "job_id": job_id}
    
    async def _check_slurm_status(self, job_id: str) -> Dict[str, Any]:
        """检查 Slurm 任务状态"""
        cmd = ["squeue", "-j", job_id, "--format=%T", "--noheader"]
        process = await asyncio.create_subprocess_exec(
            *cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE
        )
        
        stdout, stderr = await process.communicate()
        
        if process.returncode == 0:
            status = stdout.decode().strip()
            return {"status": status.lower(), "job_id": job_id}
        else:
            # 任务可能已完成
            return {"status": "completed", "job_id": job_id}
    
    async def _check_ssh_status(self, job_id: str) -> Dict[str, Any]:
        """检查 SSH 任务状态"""
        # 简化实现：检查进程是否存在
        pid = job_id.replace("ssh_", "")
        ssh_host = self.ssh_config.get("host")
        ssh_user = self.ssh_config.get("username")
        key_path = self.ssh_config.get("key_path", "~/.ssh/id_rsa")
        
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.connect(ssh_host, username=ssh_user, key_filename=os.path.expanduser(key_path))
        
        try:
            stdin, stdout, stderr = ssh.exec_command(f"ps -p {pid} > /dev/null 2>&1 && echo 'running' || echo 'completed'")
            status = stdout.read().decode().strip()
            return {"status": status, "job_id": job_id}
        finally:
            ssh.close()


def create_dispatcher_from_config(config_path: str = "config/settings.yaml") -> TaskDispatcher:
    """从配置文件创建任务分发器"""
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    dispatcher_config = config.get("dispatcher", {})
    return TaskDispatcher(dispatcher_config)


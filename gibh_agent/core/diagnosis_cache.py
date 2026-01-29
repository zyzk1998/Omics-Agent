"""
诊断结果缓存管理器
用于保存和加载数据诊断结果，避免重复诊断大型文件
"""
import json
import logging
from pathlib import Path
from typing import Dict, Any, Optional
import hashlib
import os

logger = logging.getLogger(__name__)


class DiagnosisCache:
    """
    诊断结果缓存管理器
    
    功能：
    1. 将诊断结果保存到原始数据文件同文件夹下
    2. 在规划阶段检查并加载已保存的诊断结果
    3. 支持缓存失效（基于文件修改时间）
    """
    
    def __init__(self, cache_dir: Optional[str] = None):
        """
        初始化缓存管理器
        
        Args:
            cache_dir: 缓存目录，如果为None则使用原始文件所在目录
        """
        self.cache_dir = cache_dir
    
    def _get_cache_path(self, file_path: str) -> Path:
        """
        获取缓存文件路径
        
        缓存文件命名规则：{原文件名}.diagnosis.json
        
        Args:
            file_path: 原始文件路径
            
        Returns:
            缓存文件路径
        """
        file_path_obj = Path(file_path)
        
        # 如果是目录，使用目录名
        if file_path_obj.is_dir():
            cache_name = f"{file_path_obj.name}.diagnosis.json"
            cache_path = file_path_obj / cache_name
        else:
            # 如果是文件，使用文件名（不含扩展名）+ .diagnosis.json
            cache_name = f"{file_path_obj.stem}.diagnosis.json"
            cache_path = file_path_obj.parent / cache_name
        
        return cache_path
    
    def _get_file_hash(self, file_path: str) -> str:
        """
        计算文件或目录的哈希值（用于验证缓存有效性）
        
        Args:
            file_path: 文件或目录路径
            
        Returns:
            哈希值（MD5）
        """
        file_path_obj = Path(file_path)
        
        if file_path_obj.is_dir():
            # 对于目录，使用目录中所有文件的大小和修改时间
            hasher = hashlib.md5()
            try:
                for item in sorted(file_path_obj.rglob('*')):
                    if item.is_file():
                        stat = item.stat()
                        hasher.update(f"{item.name}:{stat.st_size}:{stat.st_mtime}".encode())
            except (PermissionError, OSError) as e:
                logger.warning(f"⚠️ 无法计算目录哈希: {e}")
                return ""
            return hasher.hexdigest()
        else:
            # 对于文件，使用文件大小和修改时间
            try:
                stat = file_path_obj.stat()
                hasher = hashlib.md5()
                hasher.update(f"{stat.st_size}:{stat.st_mtime}".encode())
                return hasher.hexdigest()
            except (OSError, PermissionError) as e:
                logger.warning(f"⚠️ 无法计算文件哈希: {e}")
                return ""
    
    def save_diagnosis(
        self,
        file_path: str,
        diagnosis_result: Dict[str, Any],
        file_metadata: Optional[Dict[str, Any]] = None
    ) -> bool:
        """
        保存诊断结果到缓存文件
        
        Args:
            file_path: 原始文件路径
            diagnosis_result: 诊断结果字典
            file_metadata: 文件元数据（可选）
            
        Returns:
            True 如果保存成功，False 否则
        """
        try:
            cache_path = self._get_cache_path(file_path)
            file_hash = self._get_file_hash(file_path)
            
            # 构建缓存数据
            cache_data = {
                "file_path": str(Path(file_path).resolve()),
                "file_hash": file_hash,
                "diagnosis_result": diagnosis_result,
                "file_metadata": file_metadata,
                "cached_at": str(Path(file_path).stat().st_mtime) if Path(file_path).exists() else None
            }
            
            # 确保目录存在
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            
            # 保存到JSON文件
            with open(cache_path, 'w', encoding='utf-8') as f:
                json.dump(cache_data, f, ensure_ascii=False, indent=2)
            
            logger.info(f"✅ [DiagnosisCache] 诊断结果已保存到: {cache_path}")
            return True
            
        except Exception as e:
            logger.error(f"❌ [DiagnosisCache] 保存诊断结果失败: {e}", exc_info=True)
            return False
    
    def load_diagnosis(self, file_path: str) -> Optional[Dict[str, Any]]:
        """
        从缓存文件加载诊断结果
        
        Args:
            file_path: 原始文件路径
            
        Returns:
            诊断结果字典，如果缓存不存在或已失效则返回 None
        """
        try:
            cache_path = self._get_cache_path(file_path)
            
            if not cache_path.exists():
                logger.debug(f"ℹ️ [DiagnosisCache] 缓存文件不存在: {cache_path}")
                return None
            
            # 读取缓存数据
            with open(cache_path, 'r', encoding='utf-8') as f:
                cache_data = json.load(f)
            
            # 验证缓存有效性
            current_hash = self._get_file_hash(file_path)
            cached_hash = cache_data.get("file_hash", "")
            
            if current_hash and cached_hash and current_hash != cached_hash:
                logger.info(f"⚠️ [DiagnosisCache] 缓存已失效（文件已修改），删除旧缓存")
                try:
                    cache_path.unlink()
                except Exception as e:
                    logger.warning(f"⚠️ 删除失效缓存失败: {e}")
                return None
            
            # 验证文件路径是否匹配
            cached_file_path = cache_data.get("file_path", "")
            if cached_file_path and str(Path(file_path).resolve()) != cached_file_path:
                # 路径不匹配，但可能是同一文件的不同路径表示，继续使用
                logger.debug(f"ℹ️ [DiagnosisCache] 文件路径不匹配，但继续使用缓存")
            
            logger.info(f"✅ [DiagnosisCache] 从缓存加载诊断结果: {cache_path}")
            return cache_data.get("diagnosis_result")
            
        except json.JSONDecodeError as e:
            logger.warning(f"⚠️ [DiagnosisCache] 缓存文件格式错误，删除: {e}")
            try:
                cache_path.unlink()
            except Exception:
                pass
            return None
        except Exception as e:
            logger.error(f"❌ [DiagnosisCache] 加载诊断结果失败: {e}", exc_info=True)
            return None
    
    def has_cache(self, file_path: str) -> bool:
        """
        检查是否存在有效的缓存
        
        Args:
            file_path: 原始文件路径
            
        Returns:
            True 如果存在有效缓存，False 否则
        """
        cache_path = self._get_cache_path(file_path)
        return cache_path.exists()
    
    def clear_cache(self, file_path: str) -> bool:
        """
        清除指定文件的缓存
        
        Args:
            file_path: 原始文件路径
            
        Returns:
            True 如果清除成功，False 否则
        """
        try:
            cache_path = self._get_cache_path(file_path)
            if cache_path.exists():
                cache_path.unlink()
                logger.info(f"✅ [DiagnosisCache] 已清除缓存: {cache_path}")
                return True
            return False
        except Exception as e:
            logger.error(f"❌ [DiagnosisCache] 清除缓存失败: {e}", exc_info=True)
            return False


# 全局缓存管理器实例
_diagnosis_cache = None


def get_diagnosis_cache(cache_dir: Optional[str] = None) -> DiagnosisCache:
    """获取全局诊断缓存管理器实例"""
    global _diagnosis_cache
    if _diagnosis_cache is None:
        _diagnosis_cache = DiagnosisCache(cache_dir)
    return _diagnosis_cache

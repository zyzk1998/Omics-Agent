#!/usr/bin/env python3
"""
文件上传系统测试脚本
测试文件上传 API 的各个功能
"""
import os
import sys
import requests
import tempfile
from pathlib import Path

# 测试配置
API_BASE_URL = os.getenv("API_BASE_URL", "http://localhost:8028")
UPLOAD_ENDPOINT = f"{API_BASE_URL}/api/upload"

def create_test_file(filename: str, content: bytes = None, size_mb: float = 0.1) -> Path:
    """创建测试文件"""
    if content is None:
        # 创建指定大小的文件
        size_bytes = int(size_mb * 1024 * 1024)
        content = b"0" * size_bytes
    
    test_file = Path(tempfile.gettempdir()) / filename
    test_file.write_bytes(content)
    return test_file

def test_single_file_upload():
    """测试单文件上传"""
    print("=" * 60)
    print("测试 1: 单文件上传")
    print("=" * 60)
    
    test_file = create_test_file("test_data.csv", b"col1,col2,col3\n1,2,3\n4,5,6")
    
    try:
        with open(test_file, 'rb') as f:
            files = {'files': ('test_data.csv', f, 'text/csv')}
            response = requests.post(UPLOAD_ENDPOINT, files=files)
        
        print(f"状态码: {response.status_code}")
        if response.status_code == 200:
            data = response.json()
            print(f"✅ 上传成功")
            print(f"   文件ID: {data.get('file_id', 'N/A')}")
            print(f"   文件名: {data.get('file_name', 'N/A')}")
            print(f"   文件路径: {data.get('file_path', 'N/A')}")
            return True
        else:
            print(f"❌ 上传失败: {response.text}")
            return False
    except Exception as e:
        print(f"❌ 错误: {e}")
        return False
    finally:
        test_file.unlink(missing_ok=True)

def test_multiple_files_upload():
    """测试多文件上传"""
    print("\n" + "=" * 60)
    print("测试 2: 多文件上传")
    print("=" * 60)
    
    test_files = [
        ("file1.csv", b"col1,col2\n1,2\n3,4"),
        ("file2.tsv", b"col1\tcol2\n5\t6\n7\t8"),
    ]
    
    files = []
    temp_files = []
    
    try:
        for filename, content in test_files:
            test_file = create_test_file(filename, content)
            temp_files.append(test_file)
            files.append(('files', (filename, open(test_file, 'rb'), 'text/plain')))
        
        response = requests.post(UPLOAD_ENDPOINT, files=files)
        
        # 关闭文件
        for _, (_, file_obj, _) in files:
            file_obj.close()
        
        print(f"状态码: {response.status_code}")
        if response.status_code == 200:
            data = response.json()
            print(f"✅ 上传成功")
            print(f"   上传文件数: {data.get('count', 0)}")
            if 'files' in data:
                for file_info in data['files']:
                    print(f"   - {file_info.get('file_name', 'N/A')}")
            return True
        else:
            print(f"❌ 上传失败: {response.text}")
            return False
    except Exception as e:
        print(f"❌ 错误: {e}")
        return False
    finally:
        for f in temp_files:
            f.unlink(missing_ok=True)

def test_file_type_validation():
    """测试文件类型验证"""
    print("\n" + "=" * 60)
    print("测试 3: 文件类型验证")
    print("=" * 60)
    
    # 测试不允许的文件类型
    test_file = create_test_file("test.exe", b"malicious content")
    
    try:
        with open(test_file, 'rb') as f:
            files = {'files': ('test.exe', f, 'application/x-msdownload')}
            response = requests.post(UPLOAD_ENDPOINT, files=files)
        
        print(f"状态码: {response.status_code}")
        if response.status_code == 400:
            print(f"✅ 文件类型验证正常（拒绝了 .exe 文件）")
            print(f"   错误信息: {response.json().get('detail', 'N/A')}")
            return True
        else:
            print(f"⚠️ 文件类型验证可能有问题（应该拒绝但返回了 {response.status_code}）")
            return False
    except Exception as e:
        print(f"❌ 错误: {e}")
        return False
    finally:
        test_file.unlink(missing_ok=True)

def test_file_size_limit():
    """测试文件大小限制"""
    print("\n" + "=" * 60)
    print("测试 4: 文件大小限制")
    print("=" * 60)
    
    # 创建一个大文件（假设限制是 100MB，我们创建一个 150MB 的文件）
    print("   创建大文件（150MB）...")
    test_file = create_test_file("large_file.csv", size_mb=150)
    
    try:
        with open(test_file, 'rb') as f:
            files = {'files': ('large_file.csv', f, 'text/csv')}
            response = requests.post(UPLOAD_ENDPOINT, files=files, timeout=30)
        
        print(f"状态码: {response.status_code}")
        if response.status_code == 413:
            print(f"✅ 文件大小限制正常（拒绝了超大文件）")
            print(f"   错误信息: {response.json().get('detail', 'N/A')}")
            return True
        else:
            print(f"⚠️ 文件大小验证可能有问题（应该拒绝但返回了 {response.status_code}）")
            return False
    except requests.exceptions.Timeout:
        print("⚠️ 请求超时（可能是文件太大导致）")
        return False
    except Exception as e:
        print(f"❌ 错误: {e}")
        return False
    finally:
        test_file.unlink(missing_ok=True)

def test_10x_genomics_detection():
    """测试 10x Genomics 数据检测"""
    print("\n" + "=" * 60)
    print("测试 5: 10x Genomics 数据检测")
    print("=" * 60)
    
    # 创建 10x Genomics 文件
    test_files = [
        ("matrix.mtx", b"%%MatrixMarket matrix coordinate real general\n%metadata\n10 10 100\n"),
        ("barcodes.tsv", b"AAACCCAAGGTTCC\nAAACCCACAGGTTC\n"),
        ("features.tsv", b"ENSG00000186092\tOR4F5\tGene Expression\n"),
    ]
    
    files = []
    temp_files = []
    
    try:
        for filename, content in test_files:
            test_file = create_test_file(filename, content)
            temp_files.append(test_file)
            files.append(('files', (filename, open(test_file, 'rb'), 'text/plain')))
        
        response = requests.post(UPLOAD_ENDPOINT, files=files)
        
        # 关闭文件
        for _, (_, file_obj, _) in files:
            file_obj.close()
        
        print(f"状态码: {response.status_code}")
        if response.status_code == 200:
            data = response.json()
            print(f"✅ 上传成功")
            if 'group_dir' in data:
                print(f"   检测到 10x Genomics 数据")
                print(f"   分组目录: {data.get('group_dir', 'N/A')}")
                return True
            else:
                print(f"⚠️ 未检测到 10x Genomics 数据分组")
                return False
        else:
            print(f"❌ 上传失败: {response.text}")
            return False
    except Exception as e:
        print(f"❌ 错误: {e}")
        return False
    finally:
        for f in temp_files:
            f.unlink(missing_ok=True)

def test_filename_sanitization():
    """测试文件名清理"""
    print("\n" + "=" * 60)
    print("测试 6: 文件名清理（安全测试）")
    print("=" * 60)
    
    # 测试危险文件名
    dangerous_names = [
        "../../etc/passwd",
        "file<script>.csv",
        "file|pipe.csv",
        "file:colon.csv",
    ]
    
    results = []
    for dangerous_name in dangerous_names:
        test_file = create_test_file("safe_test.csv", b"test data")
        
        try:
            with open(test_file, 'rb') as f:
                files = {'files': (dangerous_name, f, 'text/csv')}
                response = requests.post(UPLOAD_ENDPOINT, files=files)
            
            if response.status_code == 200:
                data = response.json()
                sanitized_name = data.get('file_name', '')
                if sanitized_name != dangerous_name:
                    print(f"✅ '{dangerous_name}' -> '{sanitized_name}' (已清理)")
                    results.append(True)
                else:
                    print(f"⚠️ '{dangerous_name}' 未被清理")
                    results.append(False)
            else:
                print(f"✅ '{dangerous_name}' 被拒绝 (状态码: {response.status_code})")
                results.append(True)
        except Exception as e:
            print(f"❌ 错误: {e}")
            results.append(False)
        finally:
            test_file.unlink(missing_ok=True)
    
    return all(results)

def main():
    """运行所有测试"""
    print("=" * 60)
    print("文件上传系统测试")
    print("=" * 60)
    print(f"API 地址: {API_BASE_URL}")
    print(f"上传端点: {UPLOAD_ENDPOINT}")
    print()
    
    # 检查服务器是否可访问
    try:
        response = requests.get(API_BASE_URL, timeout=5)
        print(f"✅ 服务器可访问 (状态码: {response.status_code})")
    except requests.exceptions.ConnectionError:
        print(f"❌ 无法连接到服务器: {API_BASE_URL}")
        print("   请确保服务器正在运行")
        return 1
    except Exception as e:
        print(f"⚠️ 服务器检查失败: {e}")
    
    print()
    
    # 运行测试
    tests = [
        ("单文件上传", test_single_file_upload),
        ("多文件上传", test_multiple_files_upload),
        ("文件类型验证", test_file_type_validation),
        ("文件大小限制", test_file_size_limit),
        ("10x Genomics 检测", test_10x_genomics_detection),
        ("文件名清理", test_filename_sanitization),
    ]
    
    results = []
    for test_name, test_func in tests:
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            print(f"❌ 测试 '{test_name}' 出错: {e}")
            results.append((test_name, False))
    
    # 总结
    print("\n" + "=" * 60)
    print("测试总结")
    print("=" * 60)
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for test_name, result in results:
        status = "✅ 通过" if result else "❌ 失败"
        print(f"{status}: {test_name}")
    
    print(f"\n总计: {passed}/{total} 测试通过")
    
    return 0 if passed == total else 1

if __name__ == "__main__":
    sys.exit(main())


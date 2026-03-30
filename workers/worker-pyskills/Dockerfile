FROM continuumio/miniconda3:latest

WORKDIR /app

# 无头环境：避免 PyMOL 在无显示设备时初始化失败
ENV QT_QPA_PLATFORM=offscreen \
    MPLBACKEND=Agg

# 配置国内镜像源加速（清华源）
RUN conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/ && \
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/ && \
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/ && \
    conda config --set show_channel_urls yes

# 一键安装所有生信依赖和 FastAPI 运行环境（conda-forge / bioconda 解析依赖）
RUN conda install -y \
    python=3.10 \
    viennarna \
    pymol-open-source \
    rdkit \
    gseapy \
    matplotlib \
    numpy \
    fastapi \
    uvicorn \
    pydantic && \
    conda clean -afy

# 先拷贝代码，再校验（SA_Score 路径需见 sa_score_compat.py）
COPY . /app/

# 构建期校验：conda 版 rdkit 可能无 rdkit.Chem.SA_Score 包，由 compat 解析
RUN python -c "import RNA; import pymol; from rdkit import Chem; import gseapy; from sa_score_compat import load_calculate_score; load_calculate_score(); print('pyskills imports OK')"

EXPOSE 8000

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]

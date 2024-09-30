FROM python:3.9-slim

WORKDIR /app

COPY . /workspace

EXPOSE 8080

WORKDIR /workspace

RUN /workspace/miniconda3/envs/myenv/bin/pip install -U nupack -f /workspace/nupack-4.0.1.12/package
RUN /workspace/miniconda3/envs/myenv/bin/pip install gunicorn biopython flask wtforms flask_wtf werkzeug google-cloud-storage tqdm numpy scipy pip matplotlib pandas jupyterlab viennaRNA joblib Bio fuzzysearch

# Make RUN commands use the new environment:
#SHELL ["conda", "run", "-n", "myenv", "/bin/bash", "-c"]

CMD ["/workspace/miniconda3/envs/myenv/bin/gunicorn", "-b", "0.0.0.0:8080", "--pythonpath", "/workspace/miniconda3/envs/myenv/bin", "tool.server:app"]

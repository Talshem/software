FROM python:3.9-slim

WORKDIR /app

COPY . /workspace

EXPOSE 8080

WORKDIR /workspace

CMD ["/workspace/miniconda3/envs/myenv/bin/gunicorn", "-b", "0.0.0.0:8080", "tool.server:app"]

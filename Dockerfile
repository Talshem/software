FROM python:3.9-slim

WORKDIR /workspace

EXPOSE 8080

CMD ["source /workspace/miniconda3/bin/activate myenv && gunicorn", "-b", "0.0.0.0:8080", "tool.python:app"]

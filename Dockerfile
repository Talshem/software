FROM python:3.9-slim

WORKDIR /app

COPY . /workspace

EXPOSE 8080

CMD ["source /app/miniconda3/bin/activate myenv && gunicorn", "-b", "0.0.0.0:8080", "tool.python:app"]

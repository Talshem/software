FROM python:3.9-slim

WORKDIR /app

COPY /workspace /app/

CMD ["gunicorn", "-b", "0.0.0.0:8080", "/workspace/tool/python.py:app"]

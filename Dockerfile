FROM python:3.9-slim

WORKDIR /app

COPY workspace/ /app/

EXPOSE 8080

CMD ["gunicorn", "-b", "0.0.0.0:8080", "tool.python:app"]

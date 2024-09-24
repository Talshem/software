FROM python:3.9-slim

WORKDIR /app

COPY /workspace /app/

# Install required packages
RUN pip install --no-cache-dir -r requirements.txt

# Command to run your app
CMD ["gunicorn", "-b", "0.0.0.0:8080", "tool/python.py:app"]

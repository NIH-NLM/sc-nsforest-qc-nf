FROM python:3.11

RUN apt-get -y update && apt-get -y install cmake

COPY context/requirements-freeze.txt .
RUN pip install -r requirements-freeze.txt && pip cache purge

COPY context/ .

CMD ["python", "/main.py"]
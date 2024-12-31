FROM --platform=linux/amd64 python:3.10-slim

WORKDIR /censo
COPY src/censo /censo/src/censo
COPY src/anmr-patch-2 /censo/src/anmr-patch-2
COPY pyproject.toml /censo

RUN pip install numpy pandas matplotlib
# set python path
ENV PYTHONPATH=/censo/src
# run app as imported module and accept args from docker run
CMD ["python", "-m", "censo"]
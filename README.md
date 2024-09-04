Todo: expand this!

Build the docker image:
```sudo docker build . -t fingerprintapp```

Run the docker container, and the Genie app:
```sudo docker run -p 10008:8000 -dit fingerprintapp /bin/bash -c "julia  -e  \"using GenieFramework; Genie.loadapp(); up();\"; tail -f /dev/null"```

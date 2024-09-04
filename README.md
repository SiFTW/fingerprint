sudo docker build . -t fingerprintapp
Docker command: sudo docker run -p 10008:8000 -dit fingerprintapp /bin/bash -c "julia  -e  \"using GenieFramework; Genie.loadapp(); up();\"; tail -f /dev/null"

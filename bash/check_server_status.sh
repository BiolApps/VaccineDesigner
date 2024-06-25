#!/bin/bash

# Configuration
URL="http://155.207.86.162"
EMAIL="andigoni@auth.gr"
SUBJECT="VaccineDesigner Server Down"
BODY="The website $URL appears to be down."
CHECK_INTERVAL=300  # Check every 3600 seconds (1 hour)

# Function to check website status
check_website() {
  HTTP_STATUS=$(curl -o /dev/null -s -w "%{http_code}" $URL)
  echo "HTTP Status: $HTTP_STATUS"

  if [ "$HTTP_STATUS" -ne 200 ]; then
  curl --url 'smtp://mail.auth.gr:587' --ssl-reqd \
  --mail-from 'andigoni@auth.gr' \
  --mail-rcpt 'andigoni@auth.gr' \
  --user 'email:pass' \
  -T <(echo -e 'From: andigoni@auth.gr\nTo: andigoni@gmail.com\nSubject: VaccineDesigner is DOWN\n\nVaccineDesigner Server is not available')
fi
}

# Infinite loop
while true; do
  check_website
  sleep $CHECK_INTERVAL
done

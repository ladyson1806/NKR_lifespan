#!/usr/bin/env python3

# standard library modules
import sys, errno, re, json, ssl
from urllib import request
from urllib.error import HTTPError
from time import sleep

# BASE_URL = "https://www.ebi.ac.uk:443/interpro/api/protein/reviewed/entry/pfam/taxonomy/uniprot/10090/?page_size=200" ##PFAM MOUSE
BASE_URL = "https://www.ebi.ac.uk:443/interpro/api/protein/both/entry/pfam/taxonomy/uniprot/10181/?page_size=200" ##PFAM NAKED MOLE RAT


def output_list():
  #disable SSL verification to avoid config issues
  context = ssl._create_unverified_context()

  next = BASE_URL
  last_page = False


  #json header
  sys.stdout.write("{ \"results\": [\n")

  attempts = 0
  while next:
    try:
      req = request.Request(next, headers={"Accept": "application/json"})
      res = request.urlopen(req, context=context)
      # If the API times out due a long running query
      if res.status == 408:
        # wait just over a minute
        sleep(61)
        # then continue this loop with the same URL
        continue
      elif res.status == 204:
        #no data so leave loop
        break
      payload = json.loads(res.read().decode())
      next = payload["next"]
      attempts = 0
      if not next:
        last_page = True
    except HTTPError as e:
      if e.code == 408:
        sleep(61)
        continue
      else:
        # If there is a different HTTP error, it wil re-try 3 times before failing
        if attempts < 3:
          attempts += 1
          sleep(61)
          continue
        else:
          sys.stderr.write("LAST URL: " + next)
          raise e

    for i, item in enumerate(payload["results"]):

      sys.stdout.write(json.dumps(item))
      # for indented output replace the above line with the following
      # sys.stdout.write(json.dumps(item, indent=4))
      # for 1 record per line uncomment the following line
      # sys.stdout.write("\n")

      if last_page and i+1 == len(payload["results"]):
        sys.stdout.write("")
      else:
        sys.stdout.write(",\n")

      # Don't overload the server, give it time before asking for more
    if next:
      sleep(1)

  #json footer
  sys.stdout.write("\n] }\n")


if __name__ == "__main__":
  output_list()

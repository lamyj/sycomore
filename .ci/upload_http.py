import base64
import glob
import os
import re
import sys

token, url = [os.environ[x] for x in ("UPLOAD_TOKEN", "UPLOAD_URL")]
scheme, host, path = re.match(r"(.+?)://([^/]+)(.*)", url).groups()

if sys.version_info < (3, 0):
    import httplib
    Connection = getattr(httplib, scheme.upper()+"Connection")
else:
    import http.client
    Connection = getattr(http.client, scheme.upper()+"Connection")

try:
    auth = base64.b64encode(token)
except TypeError:
    auth = base64.b64encode(token.encode()).decode()

for wheel in glob.glob("dist/*whl"):
    print(wheel)
    data = open(wheel, "rb").read()
    connection = Connection(host)
    connection.request(
        "PUT", os.path.join(path, os.path.basename(wheel)), data,
        {"Authorization": "Basic "+auth})
    response = connection.getresponse()
    print(response.status, response.reason)
    print(response.read())
    connection.close()

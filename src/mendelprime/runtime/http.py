"""Process-scoped HTTPS compatibility; certificate verification is on by default."""

import os
import ssl
import urllib.request


def tls_context():
    if os.environ.get("MENDELPRIME_VERIFY_TLS") != "0":
        return None
    context = ssl.SSLContext(ssl.PROTOCOL_TLS_CLIENT)
    context.check_hostname = False
    context.verify_mode = ssl.CERT_NONE
    return context


def urlopen(*args, **kwargs):
    context = tls_context()
    if context is not None:
        kwargs["context"] = context
    return urllib.request.urlopen(*args, **kwargs)

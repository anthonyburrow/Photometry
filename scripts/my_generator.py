from email.generator import Generator


class MyGenerator(Generator):
    def __init__(self, fp, root=True):
        Generator.__init__(self, fp, mangle_from_=False,
                           maxheaderlen=0)
        self.root = root

    def _write_headers(self, msg):
        # We don't want to write the top-level headers;
        # they go into Request(headers) instead.
        if self.root:
            return
        # We need to use \r\n line-terminator, but Generator
        # doesn't provide the flexibility to override, so we
        # have to copy-n-paste-n-modify.
        for h, v in msg.items():
            print >> self._fp, ('%s: %s\r\n' % (h, v)),
        # A blank line always separates headers from body
        print >> self._fp, '\r\n',

    # The _write_multipart method calls "clone" for the
    # subparts.  We hijack that, setting root=False
    def clone(self, fp):
        return MyGenerator(fp, root=False)

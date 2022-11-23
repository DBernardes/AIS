# General tests of the AIS

from AIS.Header import header

from AIS.Header import Header
from tests.AIS_spectral_response_curves import ccd_operation_mode
hdr = Header(ccd_operation_mode, -70, 1)
header = hdr.create_header()
print(repr(header))

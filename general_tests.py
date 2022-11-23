from AIS.Point_Spread_Function import Point_Spread_Function
from tests.AIS_spectral_response_curves import ccd_operation_mode
import matplotlib.pyplot as plt

psf = Point_Spread_Function(ccd_operation_mode, channel=1)
image = psf.create_star_image((50, 50), 1000, 1000)
plt.imshow(image)
plt.show()

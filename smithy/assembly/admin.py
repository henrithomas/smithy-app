from django.contrib import admin
from .models import (
                GibsonAssembly, 
                GoldenGateAssembly,
                GibsonPrimer, 
                GoldenGatePrimer)

admin.site.register([GibsonAssembly, GoldenGateAssembly, GibsonPrimer, GoldenGatePrimer])


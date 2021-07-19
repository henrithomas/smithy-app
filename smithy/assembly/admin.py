from django.contrib import admin
from .models import (
                GibsonAssembly, 
                GoldenGateAssembly,
                GibsonPart,
                GoldenGatePart,
                GibsonPrimer,
                GoldenGatePrimer)

admin.site.register([
                    GibsonAssembly, 
                    GoldenGateAssembly,
                    GibsonPart,
                    GoldenGatePart,
                    GibsonPrimer,
                    GoldenGatePrimer])


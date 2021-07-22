from django.contrib import admin
from .models import (
    GibsonAssembly, 
    GoldenGateAssembly,
    GibsonPart,
    GoldenGatePart,
    GibsonPrimer,
    GoldenGatePrimer,
    GibsonSolution,
    GoldenGateSolution
)

admin.site.register([
    GibsonAssembly, 
    GoldenGateAssembly,
    GibsonPart,
    GoldenGatePart,
    GibsonPrimer,
    GoldenGatePrimer,
    GibsonSolution,
    GoldenGateSolution
])


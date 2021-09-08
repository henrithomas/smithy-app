from django.contrib import admin
from .models import (
    BioBricksAssembly,
    BioBricksPart,
    BioBricksPrimer,
    BioBricksSolution,
    GibsonAssembly, 
    GibsonPart,
    GibsonPrimer,
    GoldenGateAssembly,
    GoldenGatePart,
    GoldenGatePrimer,
    GibsonSolution,
    GoldenGateSolution,
    PCRAssembly,
    PCRSolution,
    PCRPart,
    PCRPrimer,
    SLICAssembly,
    SLICSolution,
    SLICPart,
    SLICPrimer
)

admin.site.register([
    BioBricksAssembly,
    BioBricksPart,
    BioBricksPrimer,
    BioBricksSolution,
    GibsonAssembly, 
    GibsonPart,
    GibsonPrimer,
    GoldenGateAssembly,
    GoldenGatePart,
    GoldenGatePrimer,
    GibsonSolution,
    GoldenGateSolution,
    PCRAssembly,
    PCRSolution,
    PCRPart,
    PCRPrimer,
    SLICAssembly,
    SLICSolution,
    SLICPart,
    SLICPrimer
])


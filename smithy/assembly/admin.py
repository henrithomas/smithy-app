from django.contrib import admin
from .models import GibsonAssembly, GoldenGateAssembly

admin.site.register([GibsonAssembly, GoldenGateAssembly])


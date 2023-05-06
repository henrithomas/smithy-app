from django.shortcuts import render
from itertools import chain
from ..models import (
    BioBricksAssembly,
    GibsonAssembly, 
    GoldenGateAssembly,
    PCRAssembly,
    SLICAssembly
)

def home(request):
    return render(request, 'assembly/home.html')

def about(request):
    return render(request, 'assembly/about.html', {'title': 'Assembly About'}) 

def assemblies_list(request):
    gibsons = GibsonAssembly.objects.all()
    goldengates = GoldenGateAssembly.objects.all()
    pcrs = PCRAssembly.objects.all()
    slics = SLICAssembly.objects.all()
    biobricks = BioBricksAssembly.objects.all()

    assemblies = sorted(
        chain(gibsons, goldengates, pcrs, slics, biobricks),
        key=lambda assembly: assembly.date_created, reverse=True
    )

    return render(request, 'assembly/assemblies_list.html', { 'assemblies': assemblies})
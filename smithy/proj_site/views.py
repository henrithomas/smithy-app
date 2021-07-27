from django.shortcuts import render

from time import sleep 

def home(request):
    return render(request, 'proj_site/home.html')

def about(request):
    return render(request, 'proj_site/about.html', {'title': 'About'})

def survey(request):
    return render(request, 'proj_site/survey.html', {'title': 'Survey'})

from django.shortcuts import render
import json 
import os

def home(request):
    return render(request, 'proj_site/home.html')

def about(request):
    return render(request, 'proj_site/about.html', {'title': 'About'})

def survey(request):
    return render(request, 'proj_site/survey.html', {'title': 'Survey'})

def glossary(request):
    with open('/home/hthoma/projects/smithy-app/smithy/media/glossary/smithy-glossary.json') as f:
        data = json.load(f)
        entries = data['entries']
    return render(request, 'proj_site/glossary.html', {'title': 'Glossary', 'entries': entries})
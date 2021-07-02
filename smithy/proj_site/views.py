from django.shortcuts import render

def home(request):
    return render(request, 'proj_site/home.html')

def about(request):
    return render(request, 'proj_site/about.html', {'title': 'About'})

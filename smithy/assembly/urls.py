from django.urls import path
from . import views

urlpatterns = [
    path('', views.home, name='assembly-home'),
    path('about/', views.about, name='assembly-about')
]
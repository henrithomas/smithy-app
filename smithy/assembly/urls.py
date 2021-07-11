from django.urls import path
from . import views

urlpatterns = [
    path('', views.home, name='assembly-home'),
    path('about/', views.about, name='assembly-about'),
    path('gibson/', views.gibson_form, name='gibson-form'),
    path('goldengate/', views.goldengate_form, name='goldengate-form')
]
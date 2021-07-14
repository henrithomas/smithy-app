from collections import namedtuple
from django.urls import path
from . import views
from .views import (
    GibsonDetailView,
    GibsonCreateView,
    GoldenGateDetailView,
    GoldenGateCreateView
)

urlpatterns = [
    path('', views.home, name='assembly-home'),
    path('about/', views.about, name='assembly-about'),
    path('gibson/<int:pk>/', GibsonDetailView.as_view(), name='gibson-detail'),
    path('gibson/new/', GibsonCreateView.as_view(), name='gibson-create'),
    path('goldengate/<int:pk>/', GoldenGateDetailView.as_view(), name='goldengate-detail'),
    path('goldengate/new/', GoldenGateCreateView.as_view(), name='goldengate-create'),
]
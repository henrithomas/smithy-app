from collections import namedtuple
from django.urls import path
from . import views
from .views import (
    GibsonDetailView,
    GibsonCreateView,
    GibsonPartDetailView,
    GibsonPrimerDetailView,
    GoldenGateDetailView,
    GoldenGateCreateView,
    GoldenGatePartDetailView,
    GoldenGatePrimerDetailView,
)

urlpatterns = [
    path('', views.home, name='assembly-home'),
    path('about/', views.about, name='assembly-about'),
    path('gibson/<int:pk>/', GibsonDetailView.as_view(), name='gibson-detail'),
    path('gibson/new/', GibsonCreateView.as_view(), name='gibson-create'),
    path('gibson/part/<int:pk>', GibsonPartDetailView.as_view(), name='gibson-part-detail'),
    path('gibson/primer/<int:pk>', GibsonPrimerDetailView.as_view(), name='gibson-primer-detail'),
    path('goldengate/<int:pk>/', GoldenGateDetailView.as_view(), name='goldengate-detail'),
    path('goldengate/new/', GoldenGateCreateView.as_view(), name='goldengate-create'),
    path('goldengate/part/<int:pk>', GoldenGatePartDetailView.as_view(), name='goldengate-part-detail'),
    path('goldengate/primer/<int:pk>', GoldenGatePrimerDetailView.as_view(), name='goldengate-primer-detail'),
]
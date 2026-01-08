"""
Westerlund: Westerlund ECM Panel Cointegration Test
=============================================

This package implements a functional approximation of the panel 
cointegration tests proposed by Westerlund (2007). It computes 
four statistics (Gt, Ga, Pt, Pa) based on unit-specific 
error-correction models.
"""

from .main import WesterlundTest

__all__ = ["Westerlund"]

__version__ = "0.1.0"
__author__ = "Bosco Hung"
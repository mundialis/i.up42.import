MODULE_TOPDIR = ../..

PGM = i.up42.import

include $(MODULE_TOPDIR)/include/Make/Script.make

python-requirements:
	pip install -r requirements.txt

default: python-requirements script

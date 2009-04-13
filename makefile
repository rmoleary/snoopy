# This is the default configuration
# Modifications should be made in config.mk file

MACHINE= $(shell uname -s)
HOSTNAME= $(shell hostname)
ifeq ($(MACHINE),Linux)
	DOMAINNAME= $(shell dnsdomainname)
endif

-include config.mk
include rules.mk

###############################################################
## Machine dependant config file                  #############
## This part is only used for the config rules    #############
###############################################################

# Default config file
CONFIG_FILE=src/def/config_default.mk

#Compilation variables for MACOS X (geo laptop)
ifeq ($(MACHINE),Darwin)
	CONFIG_FILE=src/def/config_mac.mk
endif

#Astro2 (DAMTP) on 32 bits
ifeq ($(HOSTNAME),astro2.damtp.cam.ac.uk)
	CONFIG_FILE=src/def/config_astro2.mk
endif

# Hyades
ifeq ($(HOSTNAME),master.hyades.private.damtp.cam.ac.uk)
	CONFIG_FILE=src/def/config_hyades.mk
endif

#Neptune and astro group coreI7
ifeq ($(HOSTNAME),neptune.damtp.cam.ac.uk)
	CONFIG_FILE=src/def/config_neptune.mk
endif

#Neptune and astro group coreI7
ifeq ($(HOSTNAME),earth.damtp.cam.ac.uk)
	CONFIG_FILE=src/def/config_earth.mk
endif

# Cambridge HPCF
ifeq ($(DOMAINNAME),moab.cluster)
	CONFIG_FILE=src/def/config_hpcf.mk
endif

#IDRIS Vargas
ifeq ($(MACHINE),AIX)
	CONFIG_FILE=src/def/config_vargas.mk
endif


#####################################################################
## Cofiguration rules (the compilation rules are in rules.mk) #######
#####################################################################
test: bench all
	@echo "******************************************"
	@echo "** Executing Benchmark *******************"
	@echo "******************************************"
	./snoopy
	@echo "******************************************"
	@echo "*** Checking results                  ****"
	@echo "******************************************"
	diff timevar src/def/timevar_bench
	@echo "******************************************"
	@echo "*** All Done  !                       ****"
	@echo "******************************************"

bench: config.mk
	@(cd src && $(MAKE) $@)
	@echo "***********************************************************"
	@echo "Initialized BENCHMARK configuration"
	@echo $(CLUSTER)
	@echo "Please type make"
	@echo "***********************************************************"



config: config.mk
	@(cd src && $(MAKE) $@)
	@echo "***********************************************************"
	@echo "Default Configuration files Initialized succesfully for"
	@echo $(CLUSTER)
	@echo "Please edit src/gvars.h and config.mk"
	@echo "***********************************************************"
	
config.mk:
	cp $(CONFIG_FILE) config.mk

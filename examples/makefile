TARGETDIR:= test

CPP:= g++
LIB_PATH:= ../Debug
INC_PATH:= ../src
INCLUDES:= -I$(TARGETDIR) -I$(INC_PATH)
CPPFLAGS:= -g3 -O0 -std=c++0x 

TEST_FILES:= $(wildcard *.mac)
TESTS:= $(patsubst %.mac, %, $(TEST_FILES))
TEST_APPS:= $(patsubst %, $(TARGETDIR)/%_app, $(TESTS))
TEST_APPS2:= $(patsubst %, $(TARGETDIR)/%_app2, $(TESTS))
TEST_APP_SRCS:= $(patsubst %, $(TARGETDIR)/%_app.cpp, $(TESTS))
TEST_APP2_SRCS:= $(patsubst %, $(TARGETDIR)/%_app2.cpp, $(TESTS))
TEST_RES:= $(patsubst %, $(TARGETDIR)/%.res, $(TESTS))
TEST_PLOTS:= $(patsubst %, $(TARGETDIR)/%_q.eps, $(TESTS))
TEST_CMPS:= $(patsubst %, $(TARGETDIR)/%_cmp.txt, $(TESTS))

all_cmps : $(TEST_CMPS)

all_apps : $(TEST_APPS2)

all : $(TEST_PLOTS)


$(TARGETDIR)/%_cmp.txt : $(TARGETDIR)/%_app2
	echo $* > $@
	cd $(TARGETDIR); ./$(<F) ../ref/$*.res | tee -a $(@F)

$(TARGETDIR)/%_q.eps : $(TARGETDIR)/%.res
	cd $(TARGETDIR)/; gnuplot $*.plt

$(TARGETDIR)/%.res : $(TARGETDIR)/%_app2
	cd $(TARGETDIR); ./$(<F)

$(TARGETDIR)/%_app : $(TARGETDIR)/%_app.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDES) $< -L$(LIB_PATH) -lCADyn -o $@

$(TARGETDIR)/%_app2 : $(TARGETDIR)/%_app2.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDES) $< -L$(LIB_PATH) -lCADyn -o $@

$(TARGETDIR)/%_app.cpp : %.mac | $(TARGETDIR)
	maxima -q --batch-string="load(\"i.lisp\")$$ load(\"../gen/cagem.mac\")$$ sysdef:prepareKinematics(\"$<\")$$ cagem2c(sysdef, \"$(TARGETDIR)\")$$ "

$(TARGETDIR)/%_app2.cpp : %.mac | $(TARGETDIR)
	maxima -q --batch-string="load(\"i.lisp\")$$ load(\"../gen/cagem.mac\")$$ sysdef:cagem(\"$<\")$$ cagem2c2(sysdef, \"$(TARGETDIR)\")$$"
	
$(TARGETDIR):
	-mkdir $(TARGETDIR)

.PHONY: clean all all_apps all_cmps
clean:
	-rm -rf $(TARGETDIR)
	
.SECONDARY: $(TEST_APPS) $(TEST_APPS2) $(TEST_APP_SRCS) $(TEST_RES) $(TEST_APP2_SRCS)

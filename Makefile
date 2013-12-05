DIR_HPP=./inc/
DIR_CPP=./src/
DIR_OBJ=./tmp/
DIR_BIN=./
DIR_MAN=./man/
DIR_TMP=./tmp/

CXXFLAGS=-O2 -Wall -pedantic -std=c++11 -I $(DIR_HPP)
LFLAGS=-lpanel -lmenu -lncurses -lboost_thread -pthread -lssl /usr/lib/libcrypto.so.1.0.0
OBJS=ChatWindow.o Contact.o ContactList.o Encoder.o InputField.o Interface.o InterfaceIndicator.o main.o MessagesList.o P2PConnection.o P2PServer.o

# search paths
vpath %.hpp $(DIR_HPP)
vpath %.cpp $(DIR_CPP)
# vpath %.o $(DIR_OBJ) # wtf - does not work

zest: $(OBJS)
	echo -e "\033[01;33m[make]\033[00;32m Linking all files..."
	echo -e "\033[01;33m[make]\033[01;36m $(addprefix $(DIR_OBJ), $^) \t\033[00;31m$(LFLAGS)\033[00m"
	$(CXX) $(LFLAGS) $(addprefix $(DIR_OBJ), $^) -o $(DIR_BIN)zest
	ln -s -f $(DIR_BIN)zest ./zest
	echo -e "\033[01;33m[make]\033[01;36m $(DIR_BIN)zest \033[00;32m has been built successfully. \033[00m"

$(OBJS): %.o: %.cpp
	@[ -d $(DIR_OBJ) ] || mkdir -p $(DIR_OBJ)
	@[ -d $(DIR_BIN) ] || mkdir -p $(DIR_BIN)
	@[ -d $(DIR_MAN) ] || mkdir -p $(DIR_MAN)
	@[ -d $(DIR_TMP) ] || mkdir -p $(DIR_TMP)
	echo -e "\033[01;33m[make]\033[01;36m $< \t\033[00;31m$(CXXFLAGS)\033[00m"
	$(CXX) $(CXXFLAGS) $< -o $(DIR_OBJ)$@ -c

doc:
	pdflatex $(DIR_DOC)src/sprawozdanie.tex
	mv $(DIR_DOC)src/*.pdf $(DIR_DOC)
	mv *.log $(DIR_TMP)
	mv *.aux $(DIR_TMP)
	mv *.pdf $(DIR_MAN)
	echo -e "\033[01;33m[make]\033[00;32m Documentation has been generated.\033[00m"

cleandoc:
	rm $(DIR_OBJ)*.aux
	rm $(DIR_OBJ)*.log
	mv $(DIR_DOC)*.pdf
	echo -e "\033[01;33m[make]\033[00;32m Documentation files have been removed.\033[00m"

clean:
	rm $(DIR_OBJ)*.o
	rm $(DIR_BIN)*
	rm ./argon
	rm ./visualisation
	echo -e "\033[01;33m[make]\033[00;32m All *.o and binary files removed.\033[00m"

.PHONY: all clean doc
.SILENT :

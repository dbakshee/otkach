PROJ = sushkov
CXX = icpc
# #XARCH=-QxHOST
MICFLAGS = -mmic
CXXFLAGS = -Wall -fopenmp -lrt -O2 -march=core-avx2
# # -opt-report -vec-report
LOGIT = 2>&1 | tee last.log ; test 0 -eq $$PIPESTATUS

# #all: sushkov.exe
# #mic: sushkov.exe

run: $(PROJ).exe
	./run.sh $(PROJ).exe a1

# OMP_NUM_THREADS=4 KMP_AFFINITY=compact,granularity=fine ./$< a1
# ./$< a1
o=o
$(PROJ).exe: sushkov-main.$o sushkov-host.$o
	$(CXX) $(CXXFLAGS) $^ -o$@

$(PROJ).exe.mic: sushkov-main.mic.$o sushkov-host.mic.$o
	$(CXX) $(CXXFLAGS) $(MICFLAGS) $^ -o$@

%.$o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o$@ $(LOGIT)

%.mic.$o: %.cpp
	$(CXX) $(CXXFLAGS) $(MICFLAGS) -c $< -o$@ $(LOGIT)

clean:
	-rm -r *.$o *.exe */*.$o */*.exe

%/.mkdir:
	mkdir -p $* && touch $@


pc: $(PROJ).pc.exe
	./$< a1

$(PROJ).pc.exe: sushkov-main.cpp sushkov-host.cpp
	cl -openmp -O2 -arch:AVX2 $^ -Fe$@

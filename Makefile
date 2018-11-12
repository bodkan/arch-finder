dir := bin

binary: $(dir) bin/arch_finder
scan:


bin/arch_finder: $(dir) src/arch_finder.cpp
	g++ -std=c++14 -Wall src/arch_finder.cpp -lboost_iostreams -o $@

clean:
	rm -rf $(dir)

$(dir):
	mkdir -p $@

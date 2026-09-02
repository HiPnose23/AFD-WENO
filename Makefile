CC = c++
CFLAGS = -O3

weno3d: weno_3d.cpp
	$(CC) $(CFLAGS) -o $@ $^

clean:
	rm weno3d

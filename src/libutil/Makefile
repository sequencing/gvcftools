#

CXXFLAGS = $(XFLAGS)

SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:%.cpp=%.o)

$(LIBUTIL): $(OBJS)
	$(AR) -csru $@ $(OBJS)

clean:
	rm -f $(LIBUTIL) $(OBJS)


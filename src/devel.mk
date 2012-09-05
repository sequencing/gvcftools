VERSION := gvcftools.hh

$(OBJS): $(VERSION)

$(VERSION): $(VERSION).in
	gitversion=$$(git describe); echo $$gitversion; sed "s/\$${VERSION}/$$gitversion/" < $@.in > $@

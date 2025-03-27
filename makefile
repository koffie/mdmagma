COMP_DIRS := $(wildcard computations/*/verify_all.sh)
.PHONY: clean copy_remote test_remote test_computations test_computations_remote $(COMP_DIRS)

clean:
	find . -type f -name "*.sig" -delete

copy_remote:
	rsync -avz --delete ./ $(ssh):/tmp/mdmagma

test:
	cd tests && magma -n v2/test_all.m

test_remote: copy_remote
	ssh $(ssh) "cd /tmp/mdmagma && make test"

test_computations_remote: copy_remote
	ssh $(ssh) "cd /tmp/mdmagma && make test_computations"

test_computations: $(COMP_DIRS)

$(COMP_DIRS):
	cd $(dir $@) && ./verify_all.sh

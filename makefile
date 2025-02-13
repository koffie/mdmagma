.PHONY: clean

clean:
	find . -type f -name "*.sig" -delete

test_remote:
	rsync -avz --delete ./ $(ssh):/tmp/mdmagma
	ssh $(ssh) "cd /tmp/mdmagma/tests && magma -n v2/test_all.m"

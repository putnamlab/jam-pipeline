
.PHONY: test


test:
	python test/test_jam_env.py
	python test/test_jam_trim.py
	python test/test_jam_bvcount.py



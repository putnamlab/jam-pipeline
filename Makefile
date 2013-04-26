
.PHONY: test all

SRC_DIR = src

project_code:
	$(MAKE) -C $(CODE_DIR)

all:
	$(MAKE) -C $(SRC_DIR) all install


test:
	python test/test_jam_env.py
	python test/test_jam_trim.py
	python test/test_jam_bvcount.py
	python test/test_jam_SNPmers.py



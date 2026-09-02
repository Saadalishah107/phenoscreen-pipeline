.PHONY: lint config test test-full clean

lint:
	python -m compileall -q bin test_data
	nextflow lint .

config:
	nextflow config .

test:
	./docker/build_images.sh
	nextflow run . -profile test,docker -resume

test-full:
	./docker/build_images.sh --full
	nextflow run . -profile test,docker -resume --dti_mode deepurpose --struct_backend mock

clean:
	rm -rf work results results_test .nextflow*

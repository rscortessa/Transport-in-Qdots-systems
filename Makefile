T=4
%.x :%.f
	gfortran -std=legacy $< -o $@
Data$(T).dat: Qdots.x
	./Qdots.x $(T)
.PHONY:
clean:
	rm *.x


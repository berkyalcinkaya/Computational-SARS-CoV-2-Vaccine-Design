with open("/Users/berk/SpikeSequence.txt", "r") as Spikeseq:
    x = Spikeseq.read()
    print(x)

base_pairs = ["G", "A", "C", "U"]
new_seq = []
for i in x:
    if i in base_pairs:
        new_seq.append(i)
print("".join(new_seq))

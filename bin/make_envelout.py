def make_envel():
    with open('envelOUT.txt') as f:
        raw_data = f.read()

    data = raw_data.split('\n \n')

    n = 0

    for i, envel in enumerate(data[5:]):
        if 'beta' in envel:
            continue

        if 'Number of critical' in envel:
            continue
        else:
            n += 1
            with open(f"envelout{n}", 'w') as w:
                w.write(envel)

if __name__ == '__main__':
    make_envel()

for i in {00..27}
do

                                 get_dynamic_contacts.py \
                                 --topology protein.pdb \
                                 --trajectory ../trajectories/rep_$i.xtc \
                                 --output contacts/cont_rep_$i.tsv \
                                 --cores 36 \
                                 --itypes all --distout
done

for i in {00..27}
do
    get_contact_frequencies.py   --input_files contacts/cont_rep_$i.tsv --output_file freqs/freqs_rep_$i.tsv
done

    def _predict_strand(self, nodes):
        total_length = sum((n[1]-n[0]) for n in nodes)
        strand_expr = {GTF.POS_STRAND: 0.0, GTF.NEG_STRAND: 0.0}
        strand_length = {GTF.POS_STRAND: 0, GTF.NEG_STRAND: 0}

        for n in nodes:
            nd = self.node_data[n]
            length = n[1] - n[0]
            frac_length = float(length) / total_length
            strand_expr['+'] += nd.exprs['+'] * frac_length
            strand_expr['-'] += nd.exprs['-'] * frac_length
            if nd.strands['+']:
                strand_length['+'] += frac_length
            if nd.strands['-']:
                strand_length['-'] += frac_length

        # if transfrag supported by stranded coverage choose strand with
        # greatest length-normalized expression level
        total_expr = sum(strand_expr.values())
        if total_expr > 0:
            if strand_expr['+'] > strand_expr['-']:
                return '+'
            elif strand_expr['-'] > strand_expr['+']:
                return '-'
        # if transfrag not supported by stranded coverage choose strand
        # with greatest support from reference transcripts
        total_strand_length = sum(strand_length.values())
        if total_strand_length > 0:
            if strand_length['+'] > strand_length['-']:
                return '+'
            elif strand_length['-'] > strand_length['+']:
                return '-'
        return '.'

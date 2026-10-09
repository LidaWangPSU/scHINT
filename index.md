# scHINT

Gusev Lab · Dana-Farber Cancer Institute

Single-cell gene expression **H**eritability with gene-by-cell-context
**INT**eractions

Estimate how much of a gene’s expression is genetic, and how much of
that genetic effect changes with cell state, cell type or perturbation,
directly from single-cell data.

[Get started](https://LidaWangPSU.github.io/scHINT/articles/scHINT.md)
[The model](https://LidaWangPSU.github.io/scHINT/articles/model.md)
[GitHub](https://github.com/LidaWangPSU/scHINT)

![](data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdib3g9IjAgMCA5ODAgMzMwIiByb2xlPSJpbWciIGFyaWEtbGFiZWw9InNjSElOVCBvdmVydmlldzogY2VsbHMgbmVzdGVkIGluIGRvbm9ycywgcGFpcnMgb2YgY2VsbHMsIEhhc2VtYW4tRWxzdG9uIHJlZ3Jlc3Npb24sIHZhcmlhbmNlIGNvbXBvbmVudHMiIGZvbnQtZmFtaWx5PSJOdW5pdG8sICYjMzk7SGVsdmV0aWNhIE5ldWUmIzM5OywgQXJpYWwsIHNhbnMtc2VyaWYiPjxzdHlsZT4udHtmb250LXNpemU6MTNweDtmaWxsOiMwMzE2MzQ7Zm9udC13ZWlnaHQ6ODAwO2xldHRlci1zcGFjaW5nOi42cHg7dGV4dC10cmFuc2Zvcm06dXBwZXJjYXNlfS5ze2ZvbnQtc2l6ZToxMnB4O2ZpbGw6IzUxNjA3YX0ubXtmb250LWZhbWlseToiRmlyYSBNb25vIixNZW5sbyxtb25vc3BhY2U7Zm9udC1zaXplOjExLjVweDtmaWxsOiMwMzE2MzR9Lmx7Zm9udC1zaXplOjEycHg7ZmlsbDojMDMxNjM0O2ZvbnQtd2VpZ2h0OjcwMH08L3N0eWxlPgo8cmVjdCB3aWR0aD0iOTgwIiBoZWlnaHQ9IjMzMCIgcng9IjE0IiBmaWxsPSIjZjZmOWZkIiBzdHJva2U9IiNkNWRlZWEiIC8+PHRleHQgY2xhc3M9InQiIHg9IjM0IiB5PSIzOCI+RGF0YTwvdGV4dD48dGV4dCBjbGFzcz0idCIgeD0iMzAwIiB5PSIzOCI+UGFpcnMgb2YgY2VsbHM8L3RleHQ+PHRleHQgY2xhc3M9InQiIHg9IjU2MCIgeT0iMzgiPlJlZ3Jlc3Npb248L3RleHQ+PHRleHQgY2xhc3M9InQiIHg9Ijc5MCIgeT0iMzgiPk91dHB1dDwvdGV4dD48cmVjdCB4PSIzNCIgeT0iNjIiIHdpZHRoPSIyMDAiIGhlaWdodD0iNjIiIHJ4PSIxMCIgZmlsbD0iI2ZmZiIgc3Ryb2tlPSIjZDVkZWVhIiAvPjx0ZXh0IGNsYXNzPSJzIiB4PSI0NCIgeT0iNzgiPmRvbm9yIGk8L3RleHQ+PHJlY3QgeD0iNDQiIHk9IjEwNSIgd2lkdGg9IjQiIGhlaWdodD0iOSIgZmlsbD0iIzlhYTdiOCIgLz48cmVjdCB4PSI1MSIgeT0iMTA0IiB3aWR0aD0iNCIgaGVpZ2h0PSIxMCIgZmlsbD0iIzlhYTdiOCIgLz48cmVjdCB4PSI1OCIgeT0iMTA3IiB3aWR0aD0iNCIgaGVpZ2h0PSI3IiBmaWxsPSIjOWFhN2I4IiAvPjxyZWN0IHg9IjY1IiB5PSI5NyIgd2lkdGg9IjQiIGhlaWdodD0iMTciIGZpbGw9IiM5YWE3YjgiIC8+PHJlY3QgeD0iNzIiIHk9IjEwMiIgd2lkdGg9IjQiIGhlaWdodD0iMTIiIGZpbGw9IiM5YWE3YjgiIC8+PHJlY3QgeD0iNzkiIHk9IjEwMSIgd2lkdGg9IjQiIGhlaWdodD0iMTMiIGZpbGw9IiM5YWE3YjgiIC8+PHJlY3QgeD0iODYiIHk9IjEwNiIgd2lkdGg9IjQiIGhlaWdodD0iOCIgZmlsbD0iIzlhYTdiOCIgLz48Y2lyY2xlIGN4PSIxMTAuMyIgY3k9IjEwNy4yIiByPSI1LjIiIGZpbGw9IiNiYWQ5ZjYiIHN0cm9rZT0iI2ZmZiIgc3Ryb2tlLXdpZHRoPSIxIj48L2NpcmNsZT48Y2lyY2xlIGN4PSIxMTQuOCIgY3k9IjkxLjMiIHI9IjUuMiIgZmlsbD0iI2IyZDRmMyIgc3Ryb2tlPSIjZmZmIiBzdHJva2Utd2lkdGg9IjEiPjwvY2lyY2xlPjxjaXJjbGUgY3g9IjExOC41IiBjeT0iMTEwLjQiIHI9IjUuMiIgZmlsbD0iI2FjZDBmMSIgc3Ryb2tlPSIjZmZmIiBzdHJva2Utd2lkdGg9IjEiPjwvY2lyY2xlPjxjaXJjbGUgY3g9IjE0MS42IiBjeT0iODQuNyIgcj0iNS4yIiBmaWxsPSIjODdiNmU0IiBzdHJva2U9IiNmZmYiIHN0cm9rZS13aWR0aD0iMSI+PC9jaXJjbGU+PGNpcmNsZSBjeD0iMTQ5LjgiIGN5PSIxMDEuOSIgcj0iNS4yIiBmaWxsPSIjNzlhZGRmIiBzdHJva2U9IiNmZmYiIHN0cm9rZS13aWR0aD0iMSI+PC9jaXJjbGU+PGNpcmNsZSBjeD0iMTY4LjQiIGN5PSI5MS4zIiByPSI1LjIiIGZpbGw9IiM1Yjk5ZDQiIHN0cm9rZT0iI2ZmZiIgc3Ryb2tlLXdpZHRoPSIxIj48L2NpcmNsZT48Y2lyY2xlIGN4PSIxNzEuNyIgY3k9IjkxLjYiIHI9IjUuMiIgZmlsbD0iIzU1OTVkMiIgc3Ryb2tlPSIjZmZmIiBzdHJva2Utd2lkdGg9IjEiPjwvY2lyY2xlPjxjaXJjbGUgY3g9IjE5OC41IiBjeT0iODguNiIgcj0iNS4yIiBmaWxsPSIjMjk3OGMyIiBzdHJva2U9IiNmZmYiIHN0cm9rZS13aWR0aD0iMSI+PC9jaXJjbGU+PGNpcmNsZSBjeD0iMjIyLjgiIGN5PSI5Mi4xIiByPSI1LjIiIGZpbGw9IiMwMTVkYjQiIHN0cm9rZT0iI2ZmZiIgc3Ryb2tlLXdpZHRoPSIxIj48L2NpcmNsZT48cmVjdCB4PSIzNCIgeT0iMTM4IiB3aWR0aD0iMjAwIiBoZWlnaHQ9IjYyIiByeD0iMTAiIGZpbGw9IiNmZmYiIHN0cm9rZT0iI2Q1ZGVlYSIgLz48dGV4dCBjbGFzcz0icyIgeD0iNDQiIHk9IjE1NCI+ZG9ub3IgajwvdGV4dD48cmVjdCB4PSI0NCIgeT0iMTczIiB3aWR0aD0iNCIgaGVpZ2h0PSIxNyIgZmlsbD0iIzlhYTdiOCIgLz48cmVjdCB4PSI1MSIgeT0iMTc5IiB3aWR0aD0iNCIgaGVpZ2h0PSIxMSIgZmlsbD0iIzlhYTdiOCIgLz48cmVjdCB4PSI1OCIgeT0iMTgzIiB3aWR0aD0iNCIgaGVpZ2h0PSI3IiBmaWxsPSIjOWFhN2I4IiAvPjxyZWN0IHg9IjY1IiB5PSIxNzUiIHdpZHRoPSI0IiBoZWlnaHQ9IjE1IiBmaWxsPSIjOWFhN2I4IiAvPjxyZWN0IHg9IjcyIiB5PSIxNzkiIHdpZHRoPSI0IiBoZWlnaHQ9IjExIiBmaWxsPSIjOWFhN2I4IiAvPjxyZWN0IHg9Ijc5IiB5PSIxNzQiIHdpZHRoPSI0IiBoZWlnaHQ9IjE2IiBmaWxsPSIjOWFhN2I4IiAvPjxyZWN0IHg9Ijg2IiB5PSIxNzgiIHdpZHRoPSI0IiBoZWlnaHQ9IjEyIiBmaWxsPSIjOWFhN2I4IiAvPjxjaXJjbGUgY3g9IjExOC40IiBjeT0iMTc2LjAiIHI9IjUuMiIgZmlsbD0iI2FkZDBmMSIgc3Ryb2tlPSIjZmZmIiBzdHJva2Utd2lkdGg9IjEiPjwvY2lyY2xlPjxjaXJjbGUgY3g9IjEyOC42IiBjeT0iMTg0LjciIHI9IjUuMiIgZmlsbD0iIzljYzVlYiIgc3Ryb2tlPSIjZmZmIiBzdHJva2Utd2lkdGg9IjEiPjwvY2lyY2xlPjxjaXJjbGUgY3g9IjE0Mi44IiBjeT0iMTgzLjciIHI9IjUuMiIgZmlsbD0iIzg0YjVlMyIgc3Ryb2tlPSIjZmZmIiBzdHJva2Utd2lkdGg9IjEiPjwvY2lyY2xlPjxjaXJjbGUgY3g9IjE2Mi45IiBjeT0iMTc0LjIiIHI9IjUuMiIgZmlsbD0iIzY0OWZkNyIgc3Ryb2tlPSIjZmZmIiBzdHJva2Utd2lkdGg9IjEiPjwvY2lyY2xlPjxjaXJjbGUgY3g9IjE2Ni43IiBjeT0iMTcxLjYiIHI9IjUuMiIgZmlsbD0iIzVkOWJkNSIgc3Ryb2tlPSIjZmZmIiBzdHJva2Utd2lkdGg9IjEiPjwvY2lyY2xlPjxjaXJjbGUgY3g9IjE3MS41IiBjeT0iMTc2LjgiIHI9IjUuMiIgZmlsbD0iIzU1OTVkMiIgc3Ryb2tlPSIjZmZmIiBzdHJva2Utd2lkdGg9IjEiPjwvY2lyY2xlPjxjaXJjbGUgY3g9IjIwOC40IiBjeT0iMTcyLjEiIHI9IjUuMiIgZmlsbD0iIzE5NmRiZCIgc3Ryb2tlPSIjZmZmIiBzdHJva2Utd2lkdGg9IjEiPjwvY2lyY2xlPjxjaXJjbGUgY3g9IjIxMy40IiBjeT0iMTY0LjUiIHI9IjUuMiIgZmlsbD0iIzExNjdiYSIgc3Ryb2tlPSIjZmZmIiBzdHJva2Utd2lkdGg9IjEiPjwvY2lyY2xlPjxjaXJjbGUgY3g9IjIxNi40IiBjeT0iMTY4LjUiIHI9IjUuMiIgZmlsbD0iIzBjNjRiOCIgc3Ryb2tlPSIjZmZmIiBzdHJva2Utd2lkdGg9IjEiPjwvY2lyY2xlPjxyZWN0IHg9IjM0IiB5PSIyMTQiIHdpZHRoPSIyMDAiIGhlaWdodD0iNjIiIHJ4PSIxMCIgZmlsbD0iI2ZmZiIgc3Ryb2tlPSIjZDVkZWVhIiAvPjx0ZXh0IGNsYXNzPSJzIiB4PSI0NCIgeT0iMjMwIj5kb25vciBrPC90ZXh0PjxyZWN0IHg9IjQ0IiB5PSIyNDgiIHdpZHRoPSI0IiBoZWlnaHQ9IjE4IiBmaWxsPSIjOWFhN2I4IiAvPjxyZWN0IHg9IjUxIiB5PSIyNjAiIHdpZHRoPSI0IiBoZWlnaHQ9IjYiIGZpbGw9IiM5YWE3YjgiIC8+PHJlY3QgeD0iNTgiIHk9IjI1OSIgd2lkdGg9IjQiIGhlaWdodD0iNyIgZmlsbD0iIzlhYTdiOCIgLz48cmVjdCB4PSI2NSIgeT0iMjYwIiB3aWR0aD0iNCIgaGVpZ2h0PSI2IiBmaWxsPSIjOWFhN2I4IiAvPjxyZWN0IHg9IjcyIiB5PSIyNTMiIHdpZHRoPSI0IiBoZWlnaHQ9IjEzIiBmaWxsPSIjOWFhN2I4IiAvPjxyZWN0IHg9Ijc5IiB5PSIyNTAiIHdpZHRoPSI0IiBoZWlnaHQ9IjE2IiBmaWxsPSIjOWFhN2I4IiAvPjxyZWN0IHg9Ijg2IiB5PSIyNTYiIHdpZHRoPSI0IiBoZWlnaHQ9IjEwIiBmaWxsPSIjOWFhN2I4IiAvPjxjaXJjbGUgY3g9IjExNS43IiBjeT0iMjQxLjEiIHI9IjUuMiIgZmlsbD0iI2IxZDNmMyIgc3Ryb2tlPSIjZmZmIiBzdHJva2Utd2lkdGg9IjEiPjwvY2lyY2xlPjxjaXJjbGUgY3g9IjEyNC44IiBjeT0iMjQ4LjIiIHI9IjUuMiIgZmlsbD0iI2EyYzllZSIgc3Ryb2tlPSIjZmZmIiBzdHJva2Utd2lkdGg9IjEiPjwvY2lyY2xlPjxjaXJjbGUgY3g9IjE1OS4yIiBjeT0iMjUyLjUiIHI9IjUuMiIgZmlsbD0iIzZhYTNkYSIgc3Ryb2tlPSIjZmZmIiBzdHJva2Utd2lkdGg9IjEiPjwvY2lyY2xlPjxjaXJjbGUgY3g9IjE2OC4yIiBjeT0iMjUzLjgiIHI9IjUuMiIgZmlsbD0iIzViOTlkNCIgc3Ryb2tlPSIjZmZmIiBzdHJva2Utd2lkdGg9IjEiPjwvY2lyY2xlPjxjaXJjbGUgY3g9IjE4MS43IiBjeT0iMjQxLjYiIHI9IjUuMiIgZmlsbD0iIzQ1OGFjYyIgc3Ryb2tlPSIjZmZmIiBzdHJva2Utd2lkdGg9IjEiPjwvY2lyY2xlPjxjaXJjbGUgY3g9IjE4My4yIiBjeT0iMjQ1LjEiIHI9IjUuMiIgZmlsbD0iIzQyODljYiIgc3Ryb2tlPSIjZmZmIiBzdHJva2Utd2lkdGg9IjEiPjwvY2lyY2xlPjxjaXJjbGUgY3g9IjE4Ni4yIiBjeT0iMjU5LjUiIHI9IjUuMiIgZmlsbD0iIzNkODVjYSIgc3Ryb2tlPSIjZmZmIiBzdHJva2Utd2lkdGg9IjEiPjwvY2lyY2xlPjxjaXJjbGUgY3g9IjE4OS4zIiBjeT0iMjU1LjgiIHI9IjUuMiIgZmlsbD0iIzM4ODJjOCIgc3Ryb2tlPSIjZmZmIiBzdHJva2Utd2lkdGg9IjEiPjwvY2lyY2xlPjxjaXJjbGUgY3g9IjIxNC4wIiBjeT0iMjQzLjgiIHI9IjUuMiIgZmlsbD0iIzEwNjdiOSIgc3Ryb2tlPSIjZmZmIiBzdHJva2Utd2lkdGg9IjEiPjwvY2lyY2xlPjx0ZXh0IGNsYXNzPSJzIiB4PSIzNCIgeT0iMjk2Ij5nZW5vdHlwZSBnIChiYXJzKSDCtyBjZWxsIHN0YXRlIGMgKGxpZ2h0IHRvIGRhcmspPC90ZXh0PjxwYXRoIGQ9Ik0yNTAgMTcwIEwyODYgMTcwIiBzdHJva2U9IiMwMDc0ZDkiIHN0cm9rZS13aWR0aD0iMiIgbWFya2VyLWVuZD0idXJsKCNhcikiIC8+PGRlZnM+PG1hcmtlciBpZD0iYXIiIG1hcmtlcndpZHRoPSI5IiBtYXJrZXJoZWlnaHQ9IjkiIHJlZng9IjciIHJlZnk9IjQuNSIgb3JpZW50PSJhdXRvIj48cGF0aCBkPSJNMCAwIEw5IDQuNSBMMCA5IHoiIGZpbGw9IiMwMDc0ZDkiIC8+PC9tYXJrZXI+PC9kZWZzPjxyZWN0IHg9IjMwMCIgeT0iNjIiIHdpZHRoPSIyMCIgaGVpZ2h0PSIyMCIgcng9IjMiIGZpbGw9IiNmNmY5ZmQiIG9wYWNpdHk9IjEuMDAiIC8+PHJlY3QgeD0iMzIyIiB5PSI2MiIgd2lkdGg9IjIwIiBoZWlnaHQ9IjIwIiByeD0iMyIgZmlsbD0iIzAwOTk5OSIgb3BhY2l0eT0iMC44NSIgLz48cmVjdCB4PSIzNDQiIHk9IjYyIiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA5OTk5IiBvcGFjaXR5PSIwLjg1IiAvPjxyZWN0IHg9IjM2NiIgeT0iNjIiIHdpZHRoPSIyMCIgaGVpZ2h0PSIyMCIgcng9IjMiIGZpbGw9IiMwMDc0ZDkiIG9wYWNpdHk9IjAuMzciIC8+PHJlY3QgeD0iMzg4IiB5PSI2MiIgd2lkdGg9IjIwIiBoZWlnaHQ9IjIwIiByeD0iMyIgZmlsbD0iIzAwNzRkOSIgb3BhY2l0eT0iMC4zNyIgLz48cmVjdCB4PSI0MTAiIHk9IjYyIiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjM3IiAvPjxyZWN0IHg9IjQzMiIgeT0iNjIiIHdpZHRoPSIyMCIgaGVpZ2h0PSIyMCIgcng9IjMiIGZpbGw9IiMwMDc0ZDkiIG9wYWNpdHk9IjAuMjAiIC8+PHJlY3QgeD0iNDU0IiB5PSI2MiIgd2lkdGg9IjIwIiBoZWlnaHQ9IjIwIiByeD0iMyIgZmlsbD0iIzAwNzRkOSIgb3BhY2l0eT0iMC4yMCIgLz48cmVjdCB4PSI0NzYiIHk9IjYyIiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjIwIiAvPjxyZWN0IHg9IjMwMCIgeT0iODQiIHdpZHRoPSIyMCIgaGVpZ2h0PSIyMCIgcng9IjMiIGZpbGw9IiMwMDk5OTkiIG9wYWNpdHk9IjAuODUiIC8+PHJlY3QgeD0iMzIyIiB5PSI4NCIgd2lkdGg9IjIwIiBoZWlnaHQ9IjIwIiByeD0iMyIgZmlsbD0iI2Y2ZjlmZCIgb3BhY2l0eT0iMS4wMCIgLz48cmVjdCB4PSIzNDQiIHk9Ijg0IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA5OTk5IiBvcGFjaXR5PSIwLjg1IiAvPjxyZWN0IHg9IjM2NiIgeT0iODQiIHdpZHRoPSIyMCIgaGVpZ2h0PSIyMCIgcng9IjMiIGZpbGw9IiMwMDc0ZDkiIG9wYWNpdHk9IjAuMzciIC8+PHJlY3QgeD0iMzg4IiB5PSI4NCIgd2lkdGg9IjIwIiBoZWlnaHQ9IjIwIiByeD0iMyIgZmlsbD0iIzAwNzRkOSIgb3BhY2l0eT0iMC4zNyIgLz48cmVjdCB4PSI0MTAiIHk9Ijg0IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjM3IiAvPjxyZWN0IHg9IjQzMiIgeT0iODQiIHdpZHRoPSIyMCIgaGVpZ2h0PSIyMCIgcng9IjMiIGZpbGw9IiMwMDc0ZDkiIG9wYWNpdHk9IjAuMjAiIC8+PHJlY3QgeD0iNDU0IiB5PSI4NCIgd2lkdGg9IjIwIiBoZWlnaHQ9IjIwIiByeD0iMyIgZmlsbD0iIzAwNzRkOSIgb3BhY2l0eT0iMC4yMCIgLz48cmVjdCB4PSI0NzYiIHk9Ijg0IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjIwIiAvPjxyZWN0IHg9IjMwMCIgeT0iMTA2IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA5OTk5IiBvcGFjaXR5PSIwLjg1IiAvPjxyZWN0IHg9IjMyMiIgeT0iMTA2IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA5OTk5IiBvcGFjaXR5PSIwLjg1IiAvPjxyZWN0IHg9IjM0NCIgeT0iMTA2IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjZjZmOWZkIiBvcGFjaXR5PSIxLjAwIiAvPjxyZWN0IHg9IjM2NiIgeT0iMTA2IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjM3IiAvPjxyZWN0IHg9IjM4OCIgeT0iMTA2IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjM3IiAvPjxyZWN0IHg9IjQxMCIgeT0iMTA2IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjM3IiAvPjxyZWN0IHg9IjQzMiIgeT0iMTA2IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjIwIiAvPjxyZWN0IHg9IjQ1NCIgeT0iMTA2IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjIwIiAvPjxyZWN0IHg9IjQ3NiIgeT0iMTA2IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjIwIiAvPjxyZWN0IHg9IjMwMCIgeT0iMTI4IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjI5IiAvPjxyZWN0IHg9IjMyMiIgeT0iMTI4IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjI5IiAvPjxyZWN0IHg9IjM0NCIgeT0iMTI4IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjI5IiAvPjxyZWN0IHg9IjM2NiIgeT0iMTI4IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjZjZmOWZkIiBvcGFjaXR5PSIxLjAwIiAvPjxyZWN0IHg9IjM4OCIgeT0iMTI4IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA5OTk5IiBvcGFjaXR5PSIwLjg1IiAvPjxyZWN0IHg9IjQxMCIgeT0iMTI4IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA5OTk5IiBvcGFjaXR5PSIwLjg1IiAvPjxyZWN0IHg9IjQzMiIgeT0iMTI4IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjM3IiAvPjxyZWN0IHg9IjQ1NCIgeT0iMTI4IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjM3IiAvPjxyZWN0IHg9IjQ3NiIgeT0iMTI4IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjM3IiAvPjxyZWN0IHg9IjMwMCIgeT0iMTUwIiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjI5IiAvPjxyZWN0IHg9IjMyMiIgeT0iMTUwIiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjI5IiAvPjxyZWN0IHg9IjM0NCIgeT0iMTUwIiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjI5IiAvPjxyZWN0IHg9IjM2NiIgeT0iMTUwIiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA5OTk5IiBvcGFjaXR5PSIwLjg1IiAvPjxyZWN0IHg9IjM4OCIgeT0iMTUwIiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjZjZmOWZkIiBvcGFjaXR5PSIxLjAwIiAvPjxyZWN0IHg9IjQxMCIgeT0iMTUwIiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA5OTk5IiBvcGFjaXR5PSIwLjg1IiAvPjxyZWN0IHg9IjQzMiIgeT0iMTUwIiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjM3IiAvPjxyZWN0IHg9IjQ1NCIgeT0iMTUwIiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjM3IiAvPjxyZWN0IHg9IjQ3NiIgeT0iMTUwIiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjM3IiAvPjxyZWN0IHg9IjMwMCIgeT0iMTcyIiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjI5IiAvPjxyZWN0IHg9IjMyMiIgeT0iMTcyIiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjI5IiAvPjxyZWN0IHg9IjM0NCIgeT0iMTcyIiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjI5IiAvPjxyZWN0IHg9IjM2NiIgeT0iMTcyIiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA5OTk5IiBvcGFjaXR5PSIwLjg1IiAvPjxyZWN0IHg9IjM4OCIgeT0iMTcyIiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA5OTk5IiBvcGFjaXR5PSIwLjg1IiAvPjxyZWN0IHg9IjQxMCIgeT0iMTcyIiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjZjZmOWZkIiBvcGFjaXR5PSIxLjAwIiAvPjxyZWN0IHg9IjQzMiIgeT0iMTcyIiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjM3IiAvPjxyZWN0IHg9IjQ1NCIgeT0iMTcyIiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjM3IiAvPjxyZWN0IHg9IjQ3NiIgeT0iMTcyIiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjM3IiAvPjxyZWN0IHg9IjMwMCIgeT0iMTk0IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjQ2IiAvPjxyZWN0IHg9IjMyMiIgeT0iMTk0IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjQ2IiAvPjxyZWN0IHg9IjM0NCIgeT0iMTk0IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjQ2IiAvPjxyZWN0IHg9IjM2NiIgeT0iMTk0IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjI5IiAvPjxyZWN0IHg9IjM4OCIgeT0iMTk0IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjI5IiAvPjxyZWN0IHg9IjQxMCIgeT0iMTk0IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjI5IiAvPjxyZWN0IHg9IjQzMiIgeT0iMTk0IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjZjZmOWZkIiBvcGFjaXR5PSIxLjAwIiAvPjxyZWN0IHg9IjQ1NCIgeT0iMTk0IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA5OTk5IiBvcGFjaXR5PSIwLjg1IiAvPjxyZWN0IHg9IjQ3NiIgeT0iMTk0IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA5OTk5IiBvcGFjaXR5PSIwLjg1IiAvPjxyZWN0IHg9IjMwMCIgeT0iMjE2IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjQ2IiAvPjxyZWN0IHg9IjMyMiIgeT0iMjE2IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjQ2IiAvPjxyZWN0IHg9IjM0NCIgeT0iMjE2IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjQ2IiAvPjxyZWN0IHg9IjM2NiIgeT0iMjE2IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjI5IiAvPjxyZWN0IHg9IjM4OCIgeT0iMjE2IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjI5IiAvPjxyZWN0IHg9IjQxMCIgeT0iMjE2IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjI5IiAvPjxyZWN0IHg9IjQzMiIgeT0iMjE2IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA5OTk5IiBvcGFjaXR5PSIwLjg1IiAvPjxyZWN0IHg9IjQ1NCIgeT0iMjE2IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjZjZmOWZkIiBvcGFjaXR5PSIxLjAwIiAvPjxyZWN0IHg9IjQ3NiIgeT0iMjE2IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA5OTk5IiBvcGFjaXR5PSIwLjg1IiAvPjxyZWN0IHg9IjMwMCIgeT0iMjM4IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjQ2IiAvPjxyZWN0IHg9IjMyMiIgeT0iMjM4IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjQ2IiAvPjxyZWN0IHg9IjM0NCIgeT0iMjM4IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjQ2IiAvPjxyZWN0IHg9IjM2NiIgeT0iMjM4IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjI5IiAvPjxyZWN0IHg9IjM4OCIgeT0iMjM4IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjI5IiAvPjxyZWN0IHg9IjQxMCIgeT0iMjM4IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIwLjI5IiAvPjxyZWN0IHg9IjQzMiIgeT0iMjM4IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA5OTk5IiBvcGFjaXR5PSIwLjg1IiAvPjxyZWN0IHg9IjQ1NCIgeT0iMjM4IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjMDA5OTk5IiBvcGFjaXR5PSIwLjg1IiAvPjxyZWN0IHg9IjQ3NiIgeT0iMjM4IiB3aWR0aD0iMjAiIGhlaWdodD0iMjAiIHJ4PSIzIiBmaWxsPSIjZjZmOWZkIiBvcGFjaXR5PSIxLjAwIiAvPjx0ZXh0IGNsYXNzPSJzIiB4PSIzMDAiIHk9IjI4MCI+cHJvZHVjdCBvZiBleHByZXNzaW9uIGZvciBldmVyeSBwYWlyIG9mIGNlbGxzPC90ZXh0PjxyZWN0IHg9IjMwMCIgeT0iMjk0IiB3aWR0aD0iMTEiIGhlaWdodD0iMTEiIHJ4PSIyIiBmaWxsPSIjMDA5OTk5IiBvcGFjaXR5PSIuODUiIC8+PHRleHQgY2xhc3M9InMiIHg9IjMxNyIgeT0iMzA0Ij5zYW1lIGRvbm9yPC90ZXh0PjxyZWN0IHg9IjQwNCIgeT0iMjk0IiB3aWR0aD0iMTEiIGhlaWdodD0iMTEiIHJ4PSIyIiBmaWxsPSIjMDA3NGQ5IiBvcGFjaXR5PSIuNCIgLz48dGV4dCBjbGFzcz0icyIgeD0iNDIxIiB5PSIzMDQiPnJlbGF0ZWQgZG9ub3JzPC90ZXh0PjxwYXRoIGQ9Ik01MTUgMTcwIEw1NDUgMTcwIiBzdHJva2U9IiMwMDc0ZDkiIHN0cm9rZS13aWR0aD0iMiIgbWFya2VyLWVuZD0idXJsKCNhcikiIC8+PHRleHQgY2xhc3M9Im0iIHg9IjU2MCIgeT0iNzIiPnk8dHNwYW4gYmFzZWxpbmUtc2hpZnQ9InN1YiIgZm9udC1zaXplPSI5Ij5pbTwvdHNwYW4+IHk8dHNwYW4gYmFzZWxpbmUtc2hpZnQ9InN1YiIgZm9udC1zaXplPSI5Ij5qbjwvdHNwYW4+IH48L3RleHQ+PHJlY3QgeD0iNTYwIiB5PSI4MSIgd2lkdGg9IjE5NiIgaGVpZ2h0PSIzNCIgcng9IjgiIGZpbGw9IiNmZmYiIHN0cm9rZT0iI2Q1ZGVlYSIgLz48dGV4dCBjbGFzcz0ibSIgeD0iNTcyIiB5PSI5NiI+MTwvdGV4dD48dGV4dCBjbGFzcz0icyIgeD0iNTcyIiB5PSIxMDkiPmludGVyY2VwdDwvdGV4dD48cmVjdCB4PSI1NjAiIHk9IjEyMyIgd2lkdGg9IjE5NiIgaGVpZ2h0PSIzNCIgcng9IjgiIGZpbGw9IiNmZmYiIHN0cm9rZT0iI2Q1ZGVlYSIgLz48dGV4dCBjbGFzcz0ibSIgeD0iNTcyIiB5PSIxMzgiPks8dHNwYW4gYmFzZWxpbmUtc2hpZnQ9InN1YiIgZm9udC1zaXplPSI5Ij5pajwvdHNwYW4+PC90ZXh0Pjx0ZXh0IGNsYXNzPSJzIiB4PSI1NzIiIHk9IjE1MSI+Z2VuZXRpYyByZWxhdGVkbmVzcyAoR1JNKTwvdGV4dD48cmVjdCB4PSI1NjAiIHk9IjE2NSIgd2lkdGg9IjE5NiIgaGVpZ2h0PSIzNCIgcng9IjgiIGZpbGw9IiNmZmYiIHN0cm9rZT0iI2Q1ZGVlYSIgLz48dGV4dCBjbGFzcz0ibSIgeD0iNTcyIiB5PSIxODAiPjEoaSA9IGopPC90ZXh0Pjx0ZXh0IGNsYXNzPSJzIiB4PSI1NzIiIHk9IjE5MyI+c2FtZSBkb25vcjwvdGV4dD48cmVjdCB4PSI1NjAiIHk9IjIwNyIgd2lkdGg9IjE5NiIgaGVpZ2h0PSIzNCIgcng9IjgiIGZpbGw9IiNmZmYiIHN0cm9rZT0iI2Q1ZGVlYSIgLz48dGV4dCBjbGFzcz0ibSIgeD0iNTcyIiB5PSIyMjIiPmM8dHNwYW4gYmFzZWxpbmUtc2hpZnQ9InN1YiIgZm9udC1zaXplPSI5Ij5pbTwvdHNwYW4+IGM8dHNwYW4gYmFzZWxpbmUtc2hpZnQ9InN1YiIgZm9udC1zaXplPSI5Ij5qbjwvdHNwYW4+PC90ZXh0Pjx0ZXh0IGNsYXNzPSJzIiB4PSI1NzIiIHk9IjIzNSI+Y2VsbCBzdGF0ZSwgY292YXJpYXRlczwvdGV4dD48cmVjdCB4PSI1NjAiIHk9IjI0OSIgd2lkdGg9IjE5NiIgaGVpZ2h0PSIzNCIgcng9IjgiIGZpbGw9IiNlNmY0ZjQiIHN0cm9rZT0iIzAwOTk5OSIgLz48dGV4dCBjbGFzcz0ibSIgeD0iNTcyIiB5PSIyNjQiPks8dHNwYW4gYmFzZWxpbmUtc2hpZnQ9InN1YiIgZm9udC1zaXplPSI5Ij5pajwvdHNwYW4+IGM8dHNwYW4gYmFzZWxpbmUtc2hpZnQ9InN1YiIgZm9udC1zaXplPSI5Ij5pbTwvdHNwYW4+IGM8dHNwYW4gYmFzZWxpbmUtc2hpZnQ9InN1YiIgZm9udC1zaXplPSI5Ij5qbjwvdHNwYW4+PC90ZXh0Pjx0ZXh0IGNsYXNzPSJzIiB4PSI1NzIiIHk9IjI3NyI+R8OXY2VsbCBzdGF0ZTwvdGV4dD48dGV4dCBjbGFzcz0icyIgeD0iNTYwIiB5PSIzMTQiPmxlYXN0IHNxdWFyZXMgb24gZG9ub3ItbGV2ZWwgc3VtczwvdGV4dD48cGF0aCBkPSJNNzY2IDE3MCBMNzg2IDE3MCIgc3Ryb2tlPSIjMDA3NGQ5IiBzdHJva2Utd2lkdGg9IjIiIG1hcmtlci1lbmQ9InVybCgjYXIpIiAvPjxyZWN0IHg9IjgwMCIgeT0iMjIxLjYiIHdpZHRoPSI0NCIgaGVpZ2h0PSI2OC40IiBmaWxsPSIjZGJlM2VlIiAvPjx0ZXh0IGNsYXNzPSJsIiB4PSI4NTQiIHk9IjI1OS44Ij5ub2lzZTwvdGV4dD48cmVjdCB4PSI4MDAiIHk9IjIwMy40IiB3aWR0aD0iNDQiIGhlaWdodD0iMTguMiIgZmlsbD0iIzlhYTdiOCIgLz48dGV4dCBjbGFzcz0ibCIgeD0iODU0IiB5PSIyMTYuNSI+Y292YXJpYXRlPC90ZXh0PjxyZWN0IHg9IjgwMCIgeT0iMTc2LjAiIHdpZHRoPSI0NCIgaGVpZ2h0PSIyNy40IiBmaWxsPSIjMGIyYTVjIiAvPjx0ZXh0IGNsYXNzPSJsIiB4PSI4NTQiIHk9IjE5My43Ij5kb25vciBJPC90ZXh0PjxyZWN0IHg9IjgwMCIgeT0iMTI1LjgiIHdpZHRoPSI0NCIgaGVpZ2h0PSI1MC4yIiBmaWxsPSIjMDA5OTk5IiAvPjx0ZXh0IGNsYXNzPSJsIiB4PSI4NTQiIHk9IjE1NC45Ij5Hw5djZWxsIHN0YXRlPC90ZXh0PjxyZWN0IHg9IjgwMCIgeT0iNjIuMCIgd2lkdGg9IjQ0IiBoZWlnaHQ9IjYzLjgiIGZpbGw9IiMwMDc0ZDkiIC8+PHRleHQgY2xhc3M9ImwiIHg9Ijg1NCIgeT0iOTcuOSI+RzwvdGV4dD48cmVjdCB4PSI4MDAiIHk9IjYyIiB3aWR0aD0iNDQiIGhlaWdodD0iMjI4IiByeD0iNiIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjZmZmIiBzdHJva2Utd2lkdGg9IjMiIC8+PHRleHQgY2xhc3M9InMiIHg9IjgwMCIgeT0iMzEyIj5leHByZXNzaW9uIHZhcmlhbmNlPC90ZXh0Pjwvc3ZnPg==)

[schint_cell()](https://LidaWangPSU.github.io/scHINT/articles/tutorial-cell-state.md)

### G×Cell-State

Main genetic and genotype × cell-state heritability from single-cell
expression of one cell type.

[schint_sample()](https://LidaWangPSU.github.io/scHINT/articles/tutorial-cell-type.md)

### G×Cell-Type

Cell-type-shared and cell-type-specific genetic variance from donor ×
cell-type pseudobulk.

[schint_pert()](https://LidaWangPSU.github.io/scHINT/articles/tutorial-perturbation.md)

### Perturbation × context

Perturbation and perturbation × cell-state / cell-type variance in
perturb-seq screens.

## 1. Introduction

**scHINT** (single-cell gene expression **H**eritability with
gene-by-cell-context **INT**eractions) quantifies how much of the
variation in a gene’s expression is explained by genetics, and how much
of that genetic effect **depends on the cellular context**. Most
expression-QTL and heritability analyses use pseudobulk profiles, which
average over cells and so cannot see genetic effects that change along a
continuous cell state. scHINT works directly on single-cell data.

**The main model is the cell-level G×Cell-State heritability model**
([`schint_cell()`](https://LidaWangPSU.github.io/scHINT/reference/schint_cell.md)):
for each gene, the variance of single-cell expression is partitioned
into

- a **main genetic** component (G), shared by all cells,
- a **G×Cell-State** component: genetic effects that vary continuously
  along cell-state axes such as principal components or pseudotime,
- a donor (individual) component, cell-state and covariate components,
  and cell-level noise.

Estimation is a Haseman-Elston (HE) regression on pairs of cells that is
assembled from donor-level sufficient statistics, so it scales to
millions of cells, and standard errors come from a donor-block jackknife
at almost no extra cost.

The same framework extends to other contexts:

| model | function | what it estimates |
|----|----|----|
| **G×Cell-State** (cell level) | [`schint_cell()`](https://LidaWangPSU.github.io/scHINT/reference/schint_cell.md) | main genetic and genotype × cell-state heritability from single-cell data |
| **G×Cell-Type** (pseudobulk) | [`schint_sample()`](https://LidaWangPSU.github.io/scHINT/reference/schint_sample.md) | cell-type-shared and cell-type-specific genetic variance from donor × cell-type expression |
| **Perturbation × context** | [`schint_pert()`](https://LidaWangPSU.github.io/scHINT/reference/schint_pert.md) | perturbation and perturbation × cell-state / cell-type variance in perturb-seq screens |

## 2. Models

Full derivations, including the sufficient-statistics algorithm and the
jackknife, are on the [Model
page](https://LidaWangPSU.github.io/scHINT/articles/model.html).

### 2.1 G×Cell-State heritability (cell level)

For a gene, let \\y\_{im}\\ be the normalized expression in cell \\m\\
(\\m=1,\dots,M_i\\) of individual \\i\\ (\\i=1,\dots,n\\),
\\g_i\in\mathbb{R}^p\\ the standardized cis-genotypes (variants within
±500 kb of the gene), \\c\_{im}\in\mathbb{R}^K\\ the cell-state features
(e.g. the leading within-cell-type principal components) and \\o\_{im}\\
further covariates:

\\ y\_{im}=g_i^\top\beta+(g_i\otimes
c\_{im})^\top\gamma+c\_{im}^\top\alpha+o\_{im}^\top\eta+e_i+\varepsilon\_{im},
\\

\\ \beta\sim N\\\left(0,\sigma_G^2p^{-1}I_p\right),\quad \gamma_k\sim
N\\\left(0,\sigma\_{G\times C,k}^2p^{-1}I_p\right),\quad e_i\sim
N(0,\sigma_I^2),\quad \varepsilon\_{im}\sim N(0,\sigma\_\varepsilon^2).
\\

\\\sigma_G^2\\ is the contribution of main cis-genetic effects and
\\\sigma\_{G\times C,k}^2\\ that of genetic effects varying along
cell-state axis \\k\\. With \\K^G\_{ij}=g_i^\top g_j/p\\ the cis-GRM,
the expected product of two different cells is

\\ \mathbb{E}(y\_{im}y\_{jn})=\sigma_G^2K^G\_{ij}+\sum\_{k=1}^{K\_{\rm
int}}\sigma\_{G\times C,k}^2K^G\_{ij}c\_{imk}c\_{jnk} +\sigma_I^2\mathbf
1(i=j)+\sum\_{k=1}^{K}\sigma\_{C,k}^2c\_{imk}c\_{jnk}+\sum\_{l=1}^{L}\sigma\_{O,l}^2o\_{iml}o\_{jnl}.
\\

scHINT regresses \\y\_{im}y\_{jn}\\ on these kernels. The baseline model

\\ y\_{im}y\_{jn}\sim 1+K^G\_{ij}+\mathbf
1(i=j)+\sum\_{k=1}^{K}c\_{imk}c\_{jnk}+\sum\_{l=1}^{L}o\_{iml}o\_{jnl}
\\

estimates the main genetic component, and the interaction model adds
\\\sum\_{k=1}^{K\_{\rm int}}K^G\_{ij}c\_{imk}c\_{jnk}\\; the
coefficients of \\K^G\_{ij}\\ and \\K^G\_{ij}c\_{imk}c\_{jnk}\\ are
\\\hat\sigma_G^2\\ and \\\hat\sigma\_{G\times C,k}^2\\, and
\\\hat\sigma\_{G\times C}^2=\sum_k\hat\sigma\_{G\times C,k}^2\\. To make
the estimates comparable to population-level heritability, which is not
diluted by cell-level noise, heritability can be defined with the
cell-level residual removed from the denominator:

\\ h_G^2=\frac{\sigma_G^2}{\sigma_G^2+\sigma\_{G\times
C}^2+\sigma_C^2+\sigma_o^2+\sigma_I^2},\qquad h\_{G\times
C}^2=\frac{\sigma\_{G\times C}^2}{\sigma_G^2+\sigma\_{G\times
C}^2+\sigma_C^2+\sigma_o^2+\sigma_I^2}. \\

A list of GRMs (e.g. MAF bins) gives one \\\sigma_G^2\\ and one
\\\sigma\_{G\times C}^2\\ per GRM.

### 2.2 G×Cell-Type heritability (pseudobulk)

For the pseudobulk expression \\y\_{it}\\ of individual \\i\\ in cell
type \\t=1,\dots,T\\:

\\
y\_{it}=g_i^\top\beta+g_i^\top\gamma_t+\alpha_t+o\_{it}^\top\eta+e_i+\varepsilon\_{it},
\\

\\ \beta\sim N\\\left(0,\sigma_G^2p^{-1}I_p\right),\quad \gamma_t\sim
N\\\left(0,\sigma\_{G\times CT}^2p^{-1}I_p\right),\quad \alpha_t\sim
N(0,\sigma_T^2),\quad e_i\sim N(0,\sigma_I^2), \\

where \\\beta\\ is the genetic effect shared across cell types,
\\\gamma_t\\ the cell-type-specific deviation, \\\alpha_t\\ the
cell-type mean and \\e_i\\ a donor effect shared across cell types. For
two observations \\(i,s)\neq(j,t)\\

\\ \mathbb{E}(y\_{is}y\_{jt})=\sigma_G^2K^G\_{ij}+\sigma\_{G\times
CT}^2K^G\_{ij}\mathbf 1(s=t) +\sigma\_{CT}^2\mathbf
1(s=t)+\sigma_I^2\mathbf 1(i=j), \\

and the HE regression uses the kernels \\K^G\_{ij}\\, \\K^G\_{ij}\mathbf
1(s=t)\\, \\\mathbf 1(s=t)\\ and \\\mathbf 1(i=j)\\ (plus covariate
kernels). \\\sigma_G^2\\ is genetic variance shared across cell types
and \\\sigma\_{G\times CT}^2\\ is cell-type-specific genetic variance.

### 2.3 Perturbation × context heritability

For a gene, let \\p_m\\ be the perturbation of cell \\m\\ and
\\c_m\in\mathbb{R}^K\\ its standardized cell-state features:

\\
y_m=b\_{p_m}+\sum\_{k=1}^{K}c\_{mk}\alpha_k+\sum\_{k=1}^{K}c\_{mk}u\_{p_mk}+\varepsilon_m,
\\

\\ b\_{p}\sim N(0,\sigma_P^2),\quad \alpha_k\sim N(0,\sigma\_{{\rm
PC},k}^2),\quad u\_{pk}\sim N(0,\sigma\_{P\times{\rm PC},k}^2),\quad
\varepsilon_m\sim N(0,\sigma\_\varepsilon^2). \\

\\\sigma_P^2\\ is perturbation-associated variation shared across cell
states and \\\sigma\_{P\times{\rm PC},k}^2\\ the variance of the
perturbation response along axis \\k\\. For two cells \\m\neq n\\ the HE
regression is

\\ y_my_n\sim 1+\mathbf 1(p_m=p_n)+\sum\_{k=1}^{K}c\_{mk}c\_{nk}+\mathbf
1(p_m=p_n)\sum\_{k=1}^{K}c\_{mk}c\_{nk}. \\

With a categorical context (cell type) \\c\_{mk}c\_{nk}\\ is replaced by
\\\mathbf 1(t_m=t_n)\\. Because genes measured in the same cells share
their pair predictors, the design cross-products are computed once for
all genes.

## 3. Install

``` r

install.packages("remotes")
remotes::install_github("LidaWangPSU/scHINT")
```

scHINT needs only base R. (`BEDMatrix` is optional, to read PLINK files
with
[`read_plink()`](https://LidaWangPSU.github.io/scHINT/reference/read_plink.md).)
The examples below use the simulated data shipped with the package:

``` r

data(schint_example)        # genotypes + cell-level and sample-level expression
ex <- schint_example
```

## 4. Usage

### 4.1 G×Cell-State heritability: `schint_cell()`

**Input**

- `data` – one row per cell (of one cell type) with expression, donor
  ID, cell-state variables and covariates.
- `y`, `id` – names of the expression and donor-ID columns.
- `geno` – donors × cis-SNPs dosage matrix (row names = donor IDs); one
  GRM is built from all SNPs. **Or** `grm` – a pre-computed GRM, or a
  list of GRMs for a multi-GRM model.
- `context` – cell-state column(s) that interact with genetics (numeric
  = continuous, factor = categorical).
- `covariates` – other covariates (main effects only).
  `jackknife = TRUE` for standard errors.

``` r

head(ex$cell[, c("donor", "state", "pc1", "pc2", "age", "sex", "batch", "GENE_A")], 3)
#>   donor       state        pc1        pc2 age sex batch     GENE_A
#> 1  D001 -0.34413687 -0.6334822 0.02127604  63   F    B1 -1.1423329
#> 2  D001  0.04947898 -0.7615181 0.84821092  63   F    B1 -0.3653666
#> 3  D001 -0.84762422  1.2939331 0.24429519  63   F    B1  0.7018750
dim(ex$geno)
#> [1] 400 250
```

**Example**

``` r

fit <- schint_cell(ex$cell, y = "GENE_A", id = "donor", geno = ex$geno,
                   context = "state",
                   covariates = c("pc1", "pc2", "age", "sex", "batch"),
                   jackknife = TRUE, n_blocks = 50)
fit
#> scHINT model: 12000 observations, 400 donors, 1 GRM
#> 
#> Variance components (variance on the standardized-expression scale):
#>             term variance     h2      se   se_h2
#>                G  0.20500 0.4140 0.04700 0.08190
#>              GxC  0.16300 0.3290 0.03890 0.04980
#>          G_total  0.36900 0.7430 0.06680 0.07080
#>                I  0.07670 0.1550 0.03580 0.07430
#>          context  0.00704 0.0142 0.00331 0.00617
#>        covariate  0.04360 0.0879 0.01130 0.02600
#>  total_explained  0.49600     NA 0.06100      NA
#> 
#> Standard errors: jackknife over 50 donor blocks
```

**Output** – `fit$summary` (grouped components) and `fit$coefficients`
(one row per regression term): `G` main genetic variance, `GxC`
G×Cell-State variance, `G_total = G + GxC`, `I` donor, `context` and
`covariate` components; `h2` population-level heritability (residual
removed from the denominator); jackknife `se*` columns. Several GRMs and
several context variables add per-GRM and per-variable rows.

``` r

# several GRMs: e.g. cis SNPs split into three MAF bins, and two cell-state axes
grms <- lapply(split_snps(ex$geno, 3, by = "maf"), make_grm)
schint_cell(ex$cell, "GENE_A", "donor", grm = grms,
            context = c("state", "subtype"), covariates = c("pc1", "pc2"))$summary[, 1:3]
#>               term    variance    fraction
#> 1                G  0.21304077  0.21304077
#> 2                1  0.06190295  0.06190295
#> 3                2  0.05429398  0.05429398
#> 4                3  0.09684385  0.09684385
#> 5              GxC  0.15004511  0.15004511
#> 6              1:C  0.05324709  0.05324709
#> 7              2:C  0.04827433  0.04827433
#> 8              3:C  0.04852369  0.04852369
#> 9          G_total  0.36308588  0.36308588
#> 10       GxC:state  0.16876151  0.16876151
#> 11     GxC:subtype -0.01871639 -0.01871639
#> 12               I  0.09303695  0.09303695
#> 13         context  0.02630493  0.02630493
#> 14       covariate  0.02027540  0.02027540
#> 15 total_explained  0.50270316  0.50270316
```

### 4.2 G×Cell-Type heritability (pseudobulk): `schint_sample()`

**Input** – `data` with one row per donor × cell type (expression, donor
ID, cell-type column, optional covariates); the same genetic input as
above; `celltype` names the cell-type column.

``` r

head(ex$sample[, c("donor", "celltype", "age", "sex", "GENE_A")], 3)
#>   donor celltype age sex     GENE_A
#> 1  D001      CT1  63   F -1.2743896
#> 2  D001      CT2  63   F -0.4401739
#> 3  D001      CT3  63   F -2.2811747
```

**Example**

``` r

fs <- schint_sample(ex$sample, y = "GENE_A", id = "donor", celltype = "celltype",
                    geno = ex$geno, covariates = c("age", "sex"), jackknife = TRUE)
fs
#> scHINT model: 2000 observations, 400 donors, 1 GRM
#> 
#> Variance components (variance on the standardized-expression scale):
#>             term variance     h2     se  se_h2
#>                G   0.1020 0.1080 0.0237 0.0248
#>              GxC   0.0566 0.0601 0.0206 0.0204
#>          G_total   0.1580 0.1680 0.0295 0.0279
#>                I   0.0228 0.0243 0.0199 0.0211
#>          context   0.7320 0.7780 0.0216 0.0215
#>        covariate   0.0278 0.0295 0.0103 0.0106
#>  total_explained   0.9410     NA 0.0362     NA
#> 
#> Standard errors: jackknife over 400 donor blocks
```

**Output** – same structure as above: `G` is genetic variance shared
across cell types, `GxC` the cell-type-specific genetic variance, `I`
the donor component shared across cell types, `context` the cell-type
mean variance. `gxc = FALSE` fits shared genetics only;
`cat_mode = "per_level"` estimates a genetic variance for each cell
type.

### 4.3 Perturbation × context heritability: `schint_pert()`

**Input** – `data` with one row per cell, one or several expression
columns (`y`), the `perturb` column, `context` (cell state and/or cell
type), `covariates`, `control` labels to drop (e.g. `"NT"`) and
optionally `perturb_group`.

``` r

data(schint_pert_example)
head(schint_pert_example[, c("perturb", "perturb_group", "celltype", "state", "GENE_A")], 3)
#>   perturb perturb_group celltype     state     GENE_A
#> 1  PERT01        group1      CT2 1.1188078 0.07107567
#> 2  PERT01        group1      CT3 1.1119961 0.14960499
#> 3  PERT01        group1      CT1 0.4430754 0.38410583
```

**Example**

``` r

fp <- schint_pert(schint_pert_example, y = c("GENE_A", "GENE_B"), perturb = "perturb",
                  context = c("state", "celltype"), covariates = c("pc1", "pc2"),
                  control = "NT", jackknife = TRUE)
fp
#> scHINT perturbation model: 3470 cells, 60 perturbations, 2 genes
#> 
#>    gene            term variance      se
#>  GENE_A               P  0.02720 0.00866
#>  GENE_A             PxC  0.06180 0.01930
#>  GENE_A         P_total  0.08900 0.02080
#>  GENE_A       PxC:state  0.05070 0.01460
#>  GENE_A    PxC:celltype  0.01110 0.01130
#>  GENE_A         context  0.09800 0.01260
#>  GENE_A       covariate  0.04050 0.00617
#>  GENE_A total_explained  0.22700 0.02510
#>  GENE_B               P  0.04350 0.01510
#>  GENE_B             PxC  0.02450 0.01220
#>  GENE_B         P_total  0.06800 0.01880
#>  GENE_B       PxC:state -0.00277 0.00293
#>  GENE_B    PxC:celltype  0.02730 0.01260
#>  GENE_B         context  0.06640 0.01250
#>  GENE_B       covariate  0.05610 0.00832
#>  GENE_B total_explained  0.19100 0.02470
#> 
#> Standard errors: jackknife over 60 perturbation blocks
```

**Output** – long-format `fp$summary` / `fp$coefficients` with one block
per gene: `P` perturbation variance shared across states, `PxC`
perturbation × context variance (split into `PxC:state` and
`PxC:celltype`), `P_total`, and per-group rows when `perturb_group` is
used.

### More

[Users
manual](https://LidaWangPSU.github.io/scHINT/articles/users-manual.html)
· [Input formats and
outputs](https://LidaWangPSU.github.io/scHINT/articles/inputs-outputs.html)
·
[Tutorials](https://LidaWangPSU.github.io/scHINT/articles/tutorial-cell-state.html)
·
[Troubleshooting](https://LidaWangPSU.github.io/scHINT/articles/troubleshooting.html)

**Citation.** Wang L, Fadil C, Wang L, Van Dyke K, Gusev A. *Single-cell
gene expression heritability informed by gene × cell-context
interactions.* (manuscript)
